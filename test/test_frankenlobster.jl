using Test
using Oceananigans.Architectures: CPU
using Oceananigans.Grids: RectilinearGrid
using Oceananigans.Fields: ConstantField
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

using OceanBioME:
    chlorophyll, conserved_tracers, PrescribedPhotosyntheticallyActiveRadiation
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    CarbonateSystem, DissolvedParticulate, Oxygen
using OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients, NitrateAmmonia
using OceanBioME.Models.NutrientsPlanktonDetritusModels: nutrient_uptake

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

function _controlled_frankenlobster(grid)
    detritus = DissolvedParticulate(
        grid;
        dissolved_remineralisation_rate=0.0,
        particulate_remineralisation_rate=(0.0, 0.0),
        sinking_speeds=(0.0, 0.0),
    )
    return FrankenLOBSTER.construct(;
        grid,
        light_attenuation=_prescribed_light(),
        nutrients=Nutrients(; nitrogen=NitrateAmmonia(; nitrification_rate=0.0)),
        detritus,
        parameters=(
            maximum_growth_rate=(P_1=1.0, P_2=1.0),
            nitrate_half_saturation=(P_1=1.0, P_2=1.0),
            alpha=(P_1=1.0, P_2=1.0),
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
end

@testset "FrankenLOBSTER public arbitrary community" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    coupled = FrankenLOBSTER.construct(;
        grid,
        light_attenuation=_prescribed_light(),
        inorganic_carbon=CarbonateSystem(),
        oxygen=Oxygen(),
        size_structure=(
            phytoplankton=(pico=[0.5], nano=[2.0]),
            zooplankton=(micro=[8.0], meso=[20.0]),
            bacterioplankton=(heterotroph=[0.4, 0.8],),
        ),
        sinking_tracers=(nano_1=0.1,),
    )
    plankton = coupled.underlying_biogeochemistry.plankton

    @test required_biogeochemical_tracers(plankton) ==
          (:nano_1, :pico_1, :meso_1, :micro_1, :heterotroph_1, :heterotroph_2)
    @test size(plankton.runtime.parameters.palatability_matrix) == (2, 4)
    @test length(unique(plankton.runtime.parameters.palatability_matrix)) > 1
    @test plankton.runtime.parameters.assimilation_matrix == fill(0.7, 2, 4)
    @test hasproperty(plankton.runtime.sinking_velocities, :nano_1)

    volume(d) = pi / 6 * d^3
    @test plankton.runtime.parameters.bacterial_maximum_uptake_rate ≈
          [1.836 / 86400 * volume(d)^0.28 for d in (0.4, 0.8)]
    @test vec(plankton.runtime.parameters.bacterial_dom_half_saturation) ≈
          [0.04284 * volume(d)^0.65 for d in (0.4, 0.8)]

    chlorophyll_field = chlorophyll(
        plankton, (tracers=(nano_1=_cell(2.0), pico_1=_cell(1.0)),)
    )
    @test chlorophyll_field[1, 1, 1] ≈ 1.31 * 3.0

    tracers = required_biogeochemical_tracers(coupled)
    @test all(t -> t in tracers, (:DOM, :sPOM, :bPOM, :DIC, :Alk, :O₂))
    groups = conserved_tracers(coupled)
    @test groups.nitrogen.nano_1 == groups.nitrogen.heterotroph_1 == 1.0
    @test groups.carbon.nano_1 == groups.carbon.heterotroph_1 == groups.carbon.DOM == 106 / 16
end

@testset "FrankenLOBSTER coupled nitrate and DOM exchange" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    bgc = _controlled_frankenlobster(grid).underlying_biogeochemistry
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)

    growth_fields = _frankenlobster_fields(; P_1=2.0)
    growth = inv(sqrt(2.0))
    @test bgc(1, 1, 1, grid, Val(:P_1), clock, growth_fields, auxiliary_fields) ≈ growth
    @test nutrient_uptake(
        1, 1, 1, grid, Val(:NO₃), bgc.plankton, bgc, growth_fields, auxiliary_fields
    ) ≈ growth
    @test nutrient_uptake(
        1, 1, 1, grid, Val(:NH₄), bgc.plankton, bgc, growth_fields, auxiliary_fields
    ) == 0.0
    @test nutrient_uptake(
        1, 1, 1, grid, bgc.plankton, bgc, growth_fields, auxiliary_fields
    ) ≈ growth

    dom_fields = _frankenlobster_fields(; DOM=3.0, H_1=2.0)
    h = bgc(1, 1, 1, grid, Val(:H_1), clock, dom_fields, auxiliary_fields)
    dom = bgc(1, 1, 1, grid, Val(:DOM), clock, dom_fields, auxiliary_fields)
    nh4 = bgc(1, 1, 1, grid, Val(:NH₄), clock, dom_fields, auxiliary_fields)
    @test [dom, h, nh4] ≈ [-3.0, 0.75, 2.25]
    @test [
        bgc(1, 1, 1, grid, Val(:sPOM), clock, dom_fields, auxiliary_fields),
        bgc(1, 1, 1, grid, Val(:bPOM), clock, dom_fields, auxiliary_fields),
    ] == [0.0, 0.0]
end
