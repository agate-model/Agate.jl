using Test
using Oceananigans.Architectures: CPU
using Oceananigans.Grids: RectilinearGrid
using Oceananigans.Fields: ConstantField
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

using OceanBioME: chlorophyll, conserved_tracers, PrescribedPhotosyntheticallyActiveRadiation
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    CarbonateSystem, DissolvedParticulate, ExplicitCalciumCarbonate, LOBSTER, Oxygen,
    nutrient_uptake
using OceanBioME.Models.NutrientsPlanktonDetritusModels.InorganicCarbonModels:
    biological_calcium_carbonate_dissolution,
    biological_calcium_carbonate_precipitation,
    particulate_calcium_carbonate_production
using OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients, NitrateAmmonia

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER

_cell(x) = fill(x, 1, 1, 1)
_prescribed_light(x=100.0) = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(x))

function _fields(;
    NO₃=1.0, NH₄=1.0, T=20.0, DOM=0.0, sPOM=0.0, bPOM=0.0,
    DIC=2000.0, Alk=2300.0, CaCO₃=0.0, S=35.0,
    P_1=0.0, P_2=0.0, Z_1=0.0, Z_2=0.0, H_1=0.0,
)
    values = (; NO₃, NH₄, T, DOM, sPOM, bPOM, DIC, Alk, CaCO₃, S, P_1, P_2, Z_1, Z_2, H_1)
    return NamedTuple{keys(values)}(Tuple(_cell(value) for value in values))
end

function _controlled(grid; parameters=(;), kwargs...)
    plankton = FrankenLOBSTER.construct(;
        grid,
        parameters=merge((
            maximum_growth_rate=(P_1=1.0, P_2=1.0),
            nitrate_half_saturation=(P_1=1.0, P_2=1.0),
            ammonia_half_saturation=(P_1=1.0, P_2=1.0),
            nitrate_ammonia_inhibition=0.1,
            light_half_saturation=1.0,
            temperature_q10=2.0,
            reference_temperature=20.0,
            phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0),
            zooplankton_excretion_rate=(Z_1=1.0, Z_2=1.0),
            zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
            maximum_predation_rate=(Z_1=0.0, Z_2=0.0),
            bacterial_maximum_uptake_rate=(H_1=2.0,),
            bacterial_dom_half_saturation=reshape([1.0], 1, 1),
            bacterial_substrate_preference=reshape([1.0], 1, 1),
            bacterial_assimilation=reshape([0.25], 1, 1),
            bacterioplankton_mortality_rate=(H_1=0.0,),
        ), parameters),
    )
    detritus = DissolvedParticulate(
        grid; dissolved_remineralisation_rate=0.0,
        particulate_remineralisation_rate=(0.0, 0.0), sinking_speeds=(0.0, 0.0),
    )
    return LOBSTER(
        grid; plankton, nutrients=Nutrients(; nitrogen=NitrateAmmonia(; nitrification_rate=0.0)),
        light_attenuation=_prescribed_light(), detritus, kwargs...,
    )
end

@testset "FrankenLOBSTER component composition and replay" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    @test_nowarn LOBSTER(grid; plankton=FrankenLOBSTER.construct())
    parameters = (
        assimilation_matrix=fill(0.65, 2, 4),
        maximum_growth_rate=(nano_1=1e-5,),
        phytoplankton_chlorophyll_ratio=1.5,
        calcium_carbonate_rain_ratio=0.2,
    )
    plankton, recipe = FrankenLOBSTER.construct_plus_recipe(;
        grid,
        size_structure=(
            phytoplankton=(pico=[0.5], nano=[2.0]),
            zooplankton=(micro=[8.0], meso=[20.0]),
            bacterioplankton=(heterotroph=[0.4, 0.8],),
        ),
        parameters,
        sinking_tracers=(nano_1=0.1,),
        open_bottom=false,
    )
    bgc = LOBSTER(
        grid; plankton, light_attenuation=_prescribed_light(),
        inorganic_carbon=CarbonateSystem(), oxygen=Oxygen(),
    )
    decoded = Agate.Construction.decode_recipe(Agate.Construction.encode_recipe(recipe))
    replayed = FrankenLOBSTER.construct(decoded; grid)

    @test required_biogeochemical_tracers(plankton) ==
          (:nano_1, :pico_1, :meso_1, :micro_1, :heterotroph_1, :heterotroph_2)
    @test size(plankton.runtime.parameters.palatability_matrix) == (2, 4)
    @test plankton.runtime.parameters.assimilation_matrix == fill(0.65, 2, 4)
    @test hasproperty(plankton.runtime.sinking_velocities, :nano_1)
    @test :Fe ∉ required_biogeochemical_tracers(bgc)
    @test all(t -> t in required_biogeochemical_tracers(bgc), (:NO₃, :NH₄, :T, :DOM, :sPOM, :bPOM))
    @test conserved_tracers(bgc).carbon.nano_1 == 6.56
    @test chlorophyll(plankton, (tracers=(nano_1=_cell(2.0), pico_1=_cell(1.0)),))[1, 1, 1] ≈ 4.5
    @test decoded == recipe
    @test recipe.definition_version == v"0.13.0"
    @test replayed.runtime.parameters == plankton.runtime.parameters
    @test replayed.calcium_carbonate_rain_ratio == plankton.calcium_carbonate_rain_ratio == 0.2
    @test_throws ArgumentError FrankenLOBSTER.construct(sinking_tracers=(P_1=0.1,))
end

@testset "FrankenLOBSTER LOBSTER physiology and exchange" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    bgc = _controlled(grid).underlying_biogeochemistry
    aux = (PAR=_cell(1.0),)
    clock = (; time=0.0)
    tendency(tracer, fields) = bgc(1, 1, 1, grid, Val(tracer), clock, fields, aux)
    uptake(tracer, fields) = nutrient_uptake(
        1, 1, 1, grid, Val(tracer), bgc.plankton, bgc, fields, aux
    )

    light = 1 - exp(-1.0)
    nitrate = _fields(; NO₃=1.0, NH₄=0.0, P_1=2.0)
    gross = uptake(:NO₃, nitrate)
    @test gross ≈ light
    @test [tendency(t, nitrate) for t in (:NO₃, :P_1, :NH₄, :DOM)] ≈
          [-gross, 0.95 * gross, 0.0375 * gross, 0.0125 * gross]

    ammonia = _fields(; NO₃=0.0, NH₄=1.0, P_1=2.0)
    @test uptake(:NH₄, ammonia) ≈ light

    mixed = _fields(; NO₃=10.0, NH₄=10.0, P_1=2.0)
    nitrate_response = 10 / 11 * exp(-1)
    ammonia_response = 10 / 11
    expected_uptake = 2 * light * (nitrate_response + ammonia_response)
    @test uptake(:NO₃, mixed) + uptake(:NH₄, mixed) ≈ expected_uptake
    @test expected_uptake > 2 * light # LOBSTER source responses are additive, not renormalized.

    warm = _fields(; NO₃=1.0, NH₄=0.0, T=30.0, P_1=2.0)
    @test tendency(:P_1, warm) ≈ 2 * tendency(:P_1, nitrate)

    excretion = _fields(; Z_1=2.0)
    @test [tendency(t, excretion) for t in (:Z_1, :NH₄, :DOM)] ≈ [-2.0, 1.0, 1.0]

    dom = _fields(; DOM=3.0, H_1=2.0)
    @test [tendency(t, dom) for t in (:DOM, :H_1, :NH₄)] ≈ [-3.0, 0.75, 2.25]
    @test [tendency(t, dom) for t in (:sPOM, :bPOM)] == [0.0, 0.0]
end

@testset "FrankenLOBSTER P-specific calcite routing" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    aux = (PAR=_cell(1.0), Ω=_cell(1.0))
    explicit_carbon() = ExplicitCalciumCarbonate(
        grid; calcium_carbonate_dissolution_rate=0.0,
        calcium_carbonate_precipitation_rate=0.0, calcium_carbonate_sinking_speed=0.0,
    )
    flux(hook, bgc, fields) = hook(1, 1, 1, grid, bgc.plankton, bgc, fields, aux)
    fluxes(bgc, fields) = [flux(hook, bgc, fields) for hook in (
        biological_calcium_carbonate_precipitation,
        particulate_calcium_carbonate_production,
        biological_calcium_carbonate_dissolution,
    )]
    scale = 0.1 * 6.56

    growth_bgc = _controlled(grid; inorganic_carbon=explicit_carbon()).underlying_biogeochemistry
    growth_fields = _fields(; NO₃=1.0, NH₄=0.0, P_1=2.0)
    retained_growth = growth_bgc(1, 1, 1, grid, Val(:P_1), (; time=0.0), growth_fields, aux)
    @test fluxes(growth_bgc, growth_fields) ≈ [scale * retained_growth, 0.0, 0.0]

    grazing_bgc = _controlled(grid; parameters=(
        maximum_growth_rate=(P_1=0.0, P_2=0.0),
        maximum_predation_rate=(Z_1=1.0, Z_2=0.0),
        zooplankton_excretion_rate=(Z_1=0.0, Z_2=0.0),
    )).underlying_biogeochemistry
    grazing_fields = _fields(; P_1=2.0, Z_1=1.0)
    grazed_P = -grazing_bgc(1, 1, 1, grid, Val(:P_1), (; time=0.0), grazing_fields, aux)
    @test fluxes(grazing_bgc, grazing_fields) ≈ [0.0, 0.7 * scale * grazed_P, 0.3 * scale * grazed_P]

    mortality_bgc = _controlled(grid; parameters=(
        maximum_growth_rate=(P_1=0.0, P_2=0.0),
        phytoplankton_mortality_rate=(P_1=1.0, P_2=0.0),
        zooplankton_excretion_rate=(Z_1=0.0, Z_2=0.0),
        zooplankton_mortality_rate=(Z_1=1.0, Z_2=0.0),
        bacterioplankton_mortality_rate=(H_1=1.0,),
    )).underlying_biogeochemistry
    mortality_fields = _fields(; P_1=2.0)
    dead_P = -mortality_bgc(1, 1, 1, grid, Val(:P_1), (; time=0.0), mortality_fields, aux)
    @test fluxes(mortality_bgc, mortality_fields) ≈ [0.0, scale * dead_P, 0.0]
    @test fluxes(mortality_bgc, _fields(; Z_1=2.0, H_1=2.0)) == zeros(3)
end
