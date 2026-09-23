using Test
using Adapt
using Oceananigans.Architectures: CPU
using Oceananigans.Grids: RectilinearGrid
using Oceananigans.Fields: ConstantField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers

using OceanBioME:
    chlorophyll, conserved_tracers, PrescribedPhotosyntheticallyActiveRadiation
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    CarbonateSystem, DissolvedParticulate, Oxygen
using OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients, NitrateAmmonia
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    dissolved_waste, inorganic_waste, nutrient_uptake, solid_waste

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER

_prescribed_light(value=100.0) =
    PrescribedPhotosyntheticallyActiveRadiation(ConstantField(value))
_cell(value) = fill(value, 1, 1, 1)

function _frankenlobster_fields(;
    NO₃=1.0, NH₄=1.0, DOM=0.0, sPOM=0.0, bPOM=0.0,
    P_1=0.0, P_2=0.0, Z_1=0.0, Z_2=0.0, H_1=0.0,
)
    return (
        NO₃=_cell(NO₃), NH₄=_cell(NH₄), DOM=_cell(DOM),
        sPOM=_cell(sPOM), bPOM=_cell(bPOM),
        P_1=_cell(P_1), P_2=_cell(P_2), Z_1=_cell(Z_1), Z_2=_cell(Z_2), H_1=_cell(H_1),
    )
end

@testset "FrankenLOBSTER public construction" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    coupled = FrankenLOBSTER.construct(;
        grid,
        light_attenuation=_prescribed_light(),
        inorganic_carbon=CarbonateSystem(),
        oxygen=Oxygen(),
        phytoplankton_chlorophyll_ratio=2.0,
    )
    npd = coupled.underlying_biogeochemistry
    plankton = npd.plankton

    @test required_biogeochemical_tracers(coupled) == (
        :NO₃, :NH₄, :P_1, :P_2, :Z_1, :Z_2, :H_1,
        :DOM, :sPOM, :bPOM, :DIC, :Alk, :O₂,
    )
    @test required_biogeochemical_auxiliary_fields(coupled) == (:PAR,)
    @test plankton.runtime.metadata.plankton_diameters == (0.6, 1.2, 6.0, 12.0, 0.6)
    @test plankton.runtime.parameters.alpha ≈ fill(0.1953 / 86400, 2)
    @test plankton.runtime.parameters.nitrate_half_saturation[1] <
          plankton.runtime.parameters.nitrate_half_saturation[2]

    groups = conserved_tracers(coupled)
    @test groups.nitrogen.P_1 == groups.nitrogen.H_1 == 1.0
    @test groups.carbon.P_1 == groups.carbon.H_1 == groups.carbon.DOM == 106 / 16

    chlorophyll_field = chlorophyll(
        plankton, (tracers=(P_1=_cell(1.0), P_2=_cell(2.0)),)
    )
    @test chlorophyll_field[1, 1, 1] ≈ 6.0

    adapted = Adapt.adapt(identity, plankton)
    @test required_biogeochemical_tracers(adapted) ==
          required_biogeochemical_tracers(plankton)
end

@testset "FrankenLOBSTER arbitrary P/Z/H community" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    coupled = FrankenLOBSTER.construct(;
        grid,
        light_attenuation=_prescribed_light(),
        size_structure=(
            phytoplankton=(pico=[0.5], nano=[2.0]),
            zooplankton=(micro=[8.0], meso=[20.0]),
            bacterioplankton=(heterotroph=[0.8, 1.6],),
        ),
    )
    plankton = coupled.underlying_biogeochemistry.plankton

    @test required_biogeochemical_tracers(plankton) ==
          (:nano_1, :pico_1, :meso_1, :micro_1, :heterotroph_1, :heterotroph_2)
    @test size(plankton.runtime.parameters.palatability_matrix) == (2, 4)
    @test size(plankton.runtime.parameters.assimilation_matrix) == (2, 4)
end

@testset "FrankenLOBSTER nitrate-only phytoplankton growth" begin
    grid = dummy_grid(Float64)
    coupled = FrankenLOBSTER.construct(;
        grid,
        light_attenuation=_prescribed_light(),
        nutrients=Nutrients(; nitrogen=NitrateAmmonia(; nitrification_rate=0.1)),
        parameters=(
            maximum_growth_rate=(P_1=1.0, P_2=1.0),
            nitrate_half_saturation=(P_1=1.0, P_2=1.0),
            alpha=(P_1=1.0, P_2=1.0),
            phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0),
            zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
            maximum_predation_rate=(Z_1=0.0, Z_2=0.0),
        ),
    )
    bgc = coupled.underlying_biogeochemistry
    fields = _frankenlobster_fields(; P_1=2.0)
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)
    growth = inv(sqrt(2.0))

    @test bgc(1, 1, 1, grid, Val(:P_1), clock, fields, auxiliary_fields) ≈ growth
    @test nutrient_uptake(
        1, 1, 1, grid, Val(:NO₃), bgc.plankton, bgc, fields, auxiliary_fields
    ) ≈ growth
    @test nutrient_uptake(
        1, 1, 1, grid, Val(:NH₄), bgc.plankton, bgc, fields, auxiliary_fields
    ) == 0.0
    @test nutrient_uptake(
        1, 1, 1, grid, bgc.plankton, bgc, fields, auxiliary_fields
    ) ≈ growth
    @test bgc(1, 1, 1, grid, Val(:NO₃), clock, fields, auxiliary_fields) ≈ 0.1 - growth
    @test bgc(1, 1, 1, grid, Val(:NH₄), clock, fields, auxiliary_fields) ≈ -0.1
end

@testset "FrankenLOBSTER Z shares one capacity across P and H" begin
    grid = dummy_grid(Float64)
    coupled = FrankenLOBSTER.construct(;
        grid,
        light_attenuation=_prescribed_light(),
        parameters=(
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0),
            zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
            bacterioplankton_mortality_rate=(H_1=0.0,),
            bacterial_maximum_uptake_rate=(H_1=0.0,),
            maximum_predation_rate=(Z_1=1.0, Z_2=0.0),
            grazing_half_saturation=(Z_1=1.0, Z_2=1.0),
            palatability_matrix=[1.0 0.0 1.0; 0.0 0.0 0.0],
            assimilation_matrix=[0.5 0.0 0.25; 0.0 0.0 0.0],
        ),
    )
    bgc = coupled.underlying_biogeochemistry
    fields = _frankenlobster_fields(; P_1=1.0, H_1=1.0, Z_1=1.0)
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)

    p = bgc(1, 1, 1, grid, Val(:P_1), clock, fields, auxiliary_fields)
    h = bgc(1, 1, 1, grid, Val(:H_1), clock, fields, auxiliary_fields)
    z = bgc(1, 1, 1, grid, Val(:Z_1), clock, fields, auxiliary_fields)
    waste = solid_waste(1, 1, 1, grid, bgc.plankton, bgc, fields, auxiliary_fields)

    @test [p, h, z, waste] ≈ [-1 / 3, -1 / 3, 1 / 4, 5 / 12]
    @test -(p + h) ≈ 2 / 3
    @test p + h + z + waste ≈ 0.0 atol=1e-14
end

@testset "FrankenLOBSTER bacterial allometry and DOM closure" begin
    sized = FrankenLOBSTER._construct_plankton(;
        size_structure=(
            phytoplankton=(P=[0.6, 1.2],),
            zooplankton=(Z=[6.0, 12.0],),
            bacterioplankton=(H=[0.4, 0.8],),
        ),
        grid=dummy_grid(Float64),
    )
    volume(d) = pi / 6 * d^3
    @test sized.runtime.parameters.bacterial_maximum_uptake_rate ≈
          [1.836 / 86400 * volume(d)^0.28 for d in (0.4, 0.8)]
    @test vec(sized.runtime.parameters.bacterial_dom_half_saturation) ≈
          [0.04284 * volume(d)^0.65 for d in (0.4, 0.8)]

    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    detritus = DissolvedParticulate(
        grid;
        dissolved_remineralisation_rate=0.0,
        particulate_remineralisation_rate=(0.0, 0.0),
        sinking_speeds=(0.0, 0.0),
    )
    coupled = FrankenLOBSTER.construct(;
        grid,
        light_attenuation=_prescribed_light(),
        nutrients=Nutrients(; nitrogen=NitrateAmmonia(; nitrification_rate=0.0)),
        detritus,
        parameters=(
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0),
            zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
            maximum_predation_rate=(Z_1=0.0, Z_2=0.0),
            bacterial_maximum_uptake_rate=(H_1=2.0,),
            bacterial_dom_half_saturation=reshape([1.0], 1, 1),
            bacterial_substrate_preference=reshape([1.0], 1, 1),
            bacterial_assimilation=reshape([0.25], 1, 1),
            bacterioplankton_mortality_rate=(H_1=0.0,),
        ),
    )
    bgc = coupled.underlying_biogeochemistry
    fields = _frankenlobster_fields(; DOM=3.0, H_1=2.0)
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)

    h = bgc(1, 1, 1, grid, Val(:H_1), clock, fields, auxiliary_fields)
    dom = bgc(1, 1, 1, grid, Val(:DOM), clock, fields, auxiliary_fields)
    nh4 = bgc(1, 1, 1, grid, Val(:NH₄), clock, fields, auxiliary_fields)

    @test [dom, h, nh4] ≈ [-3.0, 0.75, 2.25]
    @test inorganic_waste(
        1, 1, 1, grid, bgc.plankton, bgc, fields, auxiliary_fields
    ) ≈ 2.25
    @test dissolved_waste(
        1, 1, 1, grid, bgc.plankton, bgc, fields, auxiliary_fields
    ) == 0.0
    @test dom + h + nh4 ≈ 0.0 atol=1e-14
    @test bgc(1, 1, 1, grid, Val(:sPOM), clock, fields, auxiliary_fields) == 0.0
    @test bgc(1, 1, 1, grid, Val(:bPOM), clock, fields, auxiliary_fields) == 0.0
end
