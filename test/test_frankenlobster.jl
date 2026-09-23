using Test
using Oceananigans.Architectures: CPU
using Oceananigans.Grids: RectilinearGrid
using Oceananigans.Fields: ConstantField
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

using OceanBioME:
    chlorophyll, conserved_tracers, PrescribedPhotosyntheticallyActiveRadiation
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    CarbonateSystem, DissolvedParticulate, ExplicitCalciumCarbonate, Oxygen,
    biological_calcium_carbonate_dissolution,
    biological_calcium_carbonate_precipitation,
    particulate_calcium_carbonate_production
using OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients, NitrateAmmonia, Fe
using OceanBioME.Models.NutrientsPlanktonDetritusModels: nutrient_uptake

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER

_prescribed_light(value=100.0) =
    PrescribedPhotosyntheticallyActiveRadiation(ConstantField(value))
_cell(value) = fill(value, 1, 1, 1)

function _frankenlobster_fields(;
    NO₃=1.0, NH₄=1.0, Fe=1.0, T=20.0, DOM=0.0, sPOM=0.0, bPOM=0.0,
    DIC=2000.0, Alk=2300.0, CaCO₃=0.0, S=35.0,
    P_1=0.0, P_2=0.0, Z_1=0.0, Z_2=0.0, H_1=0.0,
)
    return (
        NO₃=_cell(NO₃), NH₄=_cell(NH₄), Fe=_cell(Fe), T=_cell(T), DOM=_cell(DOM),
        sPOM=_cell(sPOM), bPOM=_cell(bPOM), DIC=_cell(DIC), Alk=_cell(Alk),
        CaCO₃=_cell(CaCO₃), S=_cell(S),
        P_1=_cell(P_1), P_2=_cell(P_2), Z_1=_cell(Z_1), Z_2=_cell(Z_2), H_1=_cell(H_1),
    )
end

function _controlled_frankenlobster(grid; parameter_overrides=(;), kwargs...)
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
        parameters=merge((
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
        ), parameter_overrides),
        kwargs...,
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
    @test groups.carbon.nano_1 == groups.carbon.heterotroph_1 == groups.carbon.DOM == 6.56
    @test plankton.carbon_ratio == 6.56
    @test plankton.calcium_carbonate_rain_ratio == 0.1
    @test plankton.zooplankton_calcium_carbonate_dissolution == 0.3
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


@testset "FrankenLOBSTER P-specific calcite routing" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    auxiliary_fields = (PAR=_cell(1.0), Ω=_cell(1.0))

    calcite_flux(hook, bgc, fields) = hook(
        1, 1, 1, grid, bgc.plankton, bgc, fields, auxiliary_fields
    )
    calcite_fluxes(bgc, fields) = [
        calcite_flux(hook, bgc, fields) for hook in (
            biological_calcium_carbonate_precipitation,
            particulate_calcium_carbonate_production,
            biological_calcium_carbonate_dissolution,
        )
    ]
    calcite_scale = 0.1 * 6.56

    growth_fields = _frankenlobster_fields(; NO₃=1.0, NH₄=0.0, Fe=1e12, P_1=2.0)
    explicit_carbon() = ExplicitCalciumCarbonate(
        grid;
        calcium_carbonate_dissolution_rate=0.0,
        calcium_carbonate_precipitation_rate=0.0,
        calcium_carbonate_sinking_speed=0.0,
    )
    explicit_bgc = _controlled_frankenlobster(
        grid; inorganic_carbon=explicit_carbon()
    ).underlying_biogeochemistry
    zero_calcite_bgc = _controlled_frankenlobster(
        grid; inorganic_carbon=explicit_carbon(), calcium_carbonate_rain_ratio=0.0
    ).underlying_biogeochemistry
    retained_growth = explicit_bgc(
        1, 1, 1, grid, Val(:P_1), (; time=0.0), growth_fields, auxiliary_fields
    )
    @test calcite_fluxes(explicit_bgc, growth_fields) ≈ [calcite_scale * retained_growth, 0.0, 0.0]

    carbonate_tendency(bgc, tracer) = bgc(
        1, 1, 1, grid, Val(tracer), (; time=0.0), growth_fields, auxiliary_fields
    )
    precipitation = calcite_flux(
        biological_calcium_carbonate_precipitation, explicit_bgc, growth_fields
    )
    @test [
        carbonate_tendency(explicit_bgc, :DIC) - carbonate_tendency(zero_calcite_bgc, :DIC),
        carbonate_tendency(explicit_bgc, :Alk) - carbonate_tendency(zero_calcite_bgc, :Alk),
        carbonate_tendency(explicit_bgc, :CaCO₃),
    ] ≈ [-precipitation, -2precipitation, 0.0]

    grazing_bgc = _controlled_frankenlobster(
        grid;
        parameter_overrides=(
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            maximum_predation_rate=(Z_1=1.0, Z_2=0.0),
            zooplankton_excretion_rate=(Z_1=0.0, Z_2=0.0),
        ),
    ).underlying_biogeochemistry
    grazing_fields = _frankenlobster_fields(; P_1=2.0, Z_1=1.0)
    grazed_P = -grazing_bgc(
        1, 1, 1, grid, Val(:P_1), (; time=0.0), grazing_fields, auxiliary_fields
    )
    @test calcite_fluxes(grazing_bgc, grazing_fields) ≈
          [0.0, calcite_scale * 0.7 * grazed_P, calcite_scale * 0.3 * grazed_P]

    mortality_bgc = _controlled_frankenlobster(
        grid;
        parameter_overrides=(
            maximum_growth_rate=(P_1=0.0, P_2=0.0),
            phytoplankton_mortality_rate=(P_1=1.0, P_2=0.0),
            zooplankton_excretion_rate=(Z_1=0.0, Z_2=0.0),
            zooplankton_mortality_rate=(Z_1=1.0, Z_2=0.0),
            bacterioplankton_mortality_rate=(H_1=1.0,),
        ),
    ).underlying_biogeochemistry
    mortality_fields = _frankenlobster_fields(; P_1=2.0)
    dead_P = -mortality_bgc(
        1, 1, 1, grid, Val(:P_1), (; time=0.0), mortality_fields, auxiliary_fields
    )
    @test calcite_fluxes(mortality_bgc, mortality_fields) ≈ [0.0, calcite_scale * dead_P, 0.0]
    @test calcite_fluxes(
        mortality_bgc, _frankenlobster_fields(; Z_1=2.0, H_1=2.0)
    ) == zeros(3)
end
