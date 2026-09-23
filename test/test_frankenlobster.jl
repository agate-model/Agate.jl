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
    Nutrients, NitrateAmmonia, Fe
using OceanBioME.Models.NutrientsPlanktonDetritusModels: nutrient_uptake

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER

_prescribed_light(value=100.0) =
    PrescribedPhotosyntheticallyActiveRadiation(ConstantField(value))
_cell(value) = fill(value, 1, 1, 1)

function _frankenlobster_fields(;
    NO₃=1.0, NH₄=1.0, Fe=1.0, T=20.0, DOM=0.0, sPOM=0.0, bPOM=0.0,
    P_1=0.0, P_2=0.0, Z_1=0.0, Z_2=0.0, H_1=0.0,
)
    return (
        NO₃=_cell(NO₃), NH₄=_cell(NH₄), Fe=_cell(Fe), T=_cell(T), DOM=_cell(DOM),
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
        nutrients=Nutrients(; nitrogen=NitrateAmmonia(; nitrification_rate=0.0), iron=Fe),
        detritus,
        parameters=(
            maximum_growth_rate=(P_1=1.0, P_2=1.0),
            nitrate_half_saturation=(P_1=1.0, P_2=1.0),
            ammonium_half_saturation=(P_1=1.0, P_2=1.0),
            iron_half_saturation=(P_1=1.0, P_2=1.0),
            ammonium_inhibition=0.1,
            temperature_q10=2.0,
            reference_temperature=20.0,
            alpha=(P_1=1.0, P_2=1.0),
            phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0),
            zooplankton_excretion_rate=(Z_1=1.0, Z_2=1.0),
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
    @test plankton.runtime.parameters.ammonium_half_saturation ≈
          0.5 .* plankton.runtime.parameters.nitrate_half_saturation
    @test plankton.runtime.parameters.iron_half_saturation == fill(2e-4, 2)
    @test plankton.runtime.parameters.temperature_q10 == 1.88
    @test plankton.runtime.parameters.reference_temperature == 20.0
    @test plankton.runtime.parameters.phytoplankton_exudation_fraction == fill(0.05, 2)
    @test plankton.runtime.parameters.ammonium_fraction_of_exudate == 0.75
    @test plankton.runtime.parameters.zooplankton_excretion_rate == fill(5.8e-7, 2)
    @test plankton.runtime.parameters.ammonium_fraction_of_zooplankton_excretion == 0.5
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
    @test all(t -> t in tracers, (:NO₃, :NH₄, :Fe, :T, :DOM, :sPOM, :bPOM, :DIC, :Alk, :O₂))
    groups = conserved_tracers(coupled)
    @test groups.nitrogen.nano_1 == groups.nitrogen.heterotroph_1 == 1.0
    @test groups.iron.nano_1 == groups.iron.heterotroph_1 == 4.6375e-5
    @test !hasproperty(groups.nitrogen, :T) && !hasproperty(groups.iron, :T)
    @test groups.carbon.nano_1 == groups.carbon.heterotroph_1 == groups.carbon.DOM == 106 / 16
end

@testset "FrankenLOBSTER coupled nutrient and DOM exchange" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    bgc = _controlled_frankenlobster(grid).underlying_biogeochemistry
    auxiliary_fields = (PAR=_cell(1.0),)
    clock = (; time=0.0)

    tendency(tracer, fields) =
        bgc(1, 1, 1, grid, Val(tracer), clock, fields, auxiliary_fields)
    uptake(tracer, fields) = nutrient_uptake(
        1, 1, 1, grid, Val(tracer), bgc.plankton, bgc, fields, auxiliary_fields
    )
    total_uptake(fields) = nutrient_uptake(
        1, 1, 1, grid, bgc.plankton, bgc, fields, auxiliary_fields
    )

    light_scale = inv(sqrt(2.0))
    nitrate_only = _frankenlobster_fields(; NO₃=1.0, NH₄=0.0, Fe=1e12, P_1=2.0)
    ammonium_only = _frankenlobster_fields(; NO₃=0.0, NH₄=1.0, Fe=1e12, P_1=2.0)
    gross_nitrate_growth = uptake(:NO₃, nitrate_only)
    @test gross_nitrate_growth ≈ light_scale
    @test tendency(:P_1, nitrate_only) ≈ 0.95 * gross_nitrate_growth
    @test tendency(:NO₃, nitrate_only) ≈ -gross_nitrate_growth
    @test tendency(:NH₄, nitrate_only) ≈ 0.0375 * gross_nitrate_growth
    @test tendency(:DOM, nitrate_only) ≈ 0.0125 * gross_nitrate_growth
    nitrate_closure = sum(
        tendency(tracer, nitrate_only) for tracer in (:NO₃, :P_1, :NH₄, :DOM)
    )
    @test isapprox(nitrate_closure, 0; atol=10eps(gross_nitrate_growth))

    @test tendency(:P_1, ammonium_only) ≈ 0.95 * light_scale
    @test uptake(:NH₄, ammonium_only) ≈ light_scale

    mixed = _frankenlobster_fields(; NO₃=10.0, NH₄=10.0, Fe=1e12, P_1=2.0)
    @test tendency(:P_1, mixed) ≈ 0.95 * sqrt(2.0)
    mixed_uptake = uptake(:NO₃, mixed) + uptake(:NH₄, mixed)
    @test mixed_uptake ≈ total_uptake(mixed)
    @test mixed_uptake ≈ sqrt(2.0)
    nitrate_without_ammonium = _frankenlobster_fields(; NO₃=10.0, NH₄=0.0, Fe=1e12, P_1=2.0)
    @test uptake(:NO₃, mixed) < uptake(:NO₃, nitrate_without_ammonium)

    iron_limited = _frankenlobster_fields(; NO₃=100.0, NH₄=0.0, Fe=1.0, P_1=2.0)
    @test tendency(:P_1, iron_limited) ≈ 0.95 * light_scale
    @test uptake(:Fe, iron_limited) ≈ light_scale * 4.6375e-5

    warm = _frankenlobster_fields(; NO₃=1.0, NH₄=0.0, Fe=1e12, T=30.0, P_1=2.0)
    @test tendency(:P_1, warm) ≈ 2 * tendency(:P_1, nitrate_only)

    excretion_fields = _frankenlobster_fields(; Z_1=2.0)
    @test [
        tendency(:Z_1, excretion_fields),
        tendency(:NH₄, excretion_fields),
        tendency(:DOM, excretion_fields),
    ] ≈ [-2.0, 1.0, 1.0]

    dom_fields = _frankenlobster_fields(; DOM=3.0, H_1=2.0)
    @test [tendency(:DOM, dom_fields), tendency(:H_1, dom_fields), tendency(:NH₄, dom_fields)] ≈
          [-3.0, 0.75, 2.25]
    @test [tendency(:sPOM, dom_fields), tendency(:bPOM, dom_fields)] == [0.0, 0.0]
end
