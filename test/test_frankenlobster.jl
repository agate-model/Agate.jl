using Test
using Oceananigans.Architectures: CPU
using Oceananigans.Grids: RectilinearGrid
using Oceananigans.Fields: ConstantField
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

using OceanBioME: chlorophyll, PrescribedPhotosyntheticallyActiveRadiation
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    DissolvedParticulate, LOBSTER, nutrient_uptake
using OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients, NitrateAmmonia

const FrankenLOBSTER = Agate.Models.FrankenLOBSTER
const _GRID = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
_cell(x) = fill(x, 1, 1, 1)
_light(x=1.0) = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(x))

function _fields(; NO₃=0.0, NH₄=0.0, T=20.0, DOM=0.0, sPOM=0.0, bPOM=0.0,
                 P_1=0.0, P_2=0.0, Z_1=0.0, Z_2=0.0, H_1=0.0)
    state = (; NO₃, NH₄, T, DOM, sPOM, bPOM, P_1, P_2, Z_1, Z_2, H_1)
    return NamedTuple{keys(state)}(map(_cell, values(state)))
end

const _CONTROLLED = (
    maximum_growth_rate=(P_1=1.0, P_2=1.0), nitrate_half_saturation=(P_1=1.0, P_2=1.0),
    ammonia_half_saturation=(P_1=1.0, P_2=1.0), nitrate_ammonia_inhibition=0.1,
    light_half_saturation=(P_1=1.0, P_2=1.0),
    temperature_q10=(P_1=2.0, P_2=2.0), reference_temperature=20.0,
    phytoplankton_mortality_rate=(P_1=0.0, P_2=0.0), maximum_predation_rate=(Z_1=0.0, Z_2=0.0),
    zooplankton_excretion_rate=(Z_1=1.0, Z_2=1.0), zooplankton_mortality_rate=(Z_1=0.0, Z_2=0.0),
    bacterial_maximum_uptake_rate=(H_1=2.0,), bacterial_dom_half_saturation=reshape([1.0], 1, 1),
    bacterial_substrate_preference=reshape([1.0], 1, 1), bacterial_assimilation=reshape([0.25], 1, 1),
    bacterioplankton_mortality_rate=(H_1=0.0,),
)

function _controlled(; parameters=(;))
    plankton = FrankenLOBSTER.construct(; grid=_GRID, parameters=merge(_CONTROLLED, parameters))
    detritus = DissolvedParticulate(
        _GRID; dissolved_remineralisation_rate=0.0,
        particulate_remineralisation_rate=(0.0, 0.0), sinking_speeds=(0.0, 0.0),
    )
    return LOBSTER(
        _GRID; plankton, nutrients=Nutrients(; nitrogen=NitrateAmmonia(; nitrification_rate=0.0)),
        light_attenuation=_light(), detritus,
    ).underlying_biogeochemistry
end

@testset "FrankenLOBSTER composition and replay" begin
    size_structure = (
        phytoplankton=(pico=[0.5], nano=[2.0]), zooplankton=(micro=[8.0], meso=[20.0]),
        bacterioplankton=(heterotroph=[0.4, 0.8],),
    )
    parameters = (
        assimilation_matrix=fill(0.65, 2, 4), maximum_growth_rate=(nano_1=1e-5,),
    )
    settings = (chlorophyll_ratio=1.5,)
    plankton, recipe = FrankenLOBSTER.construct_plus_recipe(;
        grid=_GRID, size_structure, parameters, settings,
        sinking_tracers=(nano_1=0.1,), open_bottom=false,
    )
    replayed = FrankenLOBSTER.construct(
        Agate.Construction.decode_recipe(Agate.Construction.encode_recipe(recipe)); grid=_GRID,
    )
    bgc = LOBSTER(_GRID; plankton)

    @test required_biogeochemical_tracers(plankton) ==
          (:nano_1, :pico_1, :meso_1, :micro_1, :heterotroph_1, :heterotroph_2)
    @test size(plankton.runtime.parameters.bacterial_dom_half_saturation) == (2, 1)
    @test all(t -> t in required_biogeochemical_tracers(bgc), (:NO₃, :NH₄, :DOM, :sPOM, :bPOM, :T))
    @test chlorophyll(plankton, (tracers=(nano_1=_cell(2.0), pico_1=_cell(1.0)),))[1, 1, 1] ≈ 4.5
    @test recipe.setting_overrides == settings
    @test Agate.Introspection.model_settings(plankton.runtime).chlorophyll_ratio == 1.5
    @test (replayed.runtime.parameters, replayed.traits) == (plankton.runtime.parameters, plankton.traits)
    @test_throws ArgumentError FrankenLOBSTER.construct(settings=(unknown=1.0,))
    @test_throws ArgumentError FrankenLOBSTER.construct(sinking_tracers=(P_1=0.1,))
end

@testset "FrankenLOBSTER LOBSTER exchange" begin
    bgc = _controlled()
    aux = (PAR=_cell(1.0),)
    tendency(tracer, fields) = bgc(1, 1, 1, _GRID, Val(tracer), (; time=0.0), fields, aux)
    uptake(tracer, fields) = nutrient_uptake(1, 1, 1, _GRID, Val(tracer), bgc.plankton, bgc, fields, aux)

    nitrate = _fields(; NO₃=1.0, P_1=2.0)
    gross = uptake(:NO₃, nitrate)
    @test [gross, tendency(:P_1, nitrate), tendency(:NH₄, nitrate), tendency(:DOM, nitrate)] ≈
          [1 - exp(-1), 0.95 * gross, 0.0375 * gross, 0.0125 * gross]

    mixed = _fields(; NO₃=10.0, NH₄=10.0, P_1=2.0)
    @test uptake(:NO₃, mixed) + uptake(:NH₄, mixed) ≈
          2 * (1 - exp(-1)) * (10 / 11 * exp(-1) + 10 / 11)
    @test [tendency(t, _fields(; Z_1=2.0)) for t in (:Z_1, :NH₄, :DOM)] ≈ [-2.0, 1.0, 1.0]
    @test [tendency(t, _fields(; DOM=3.0, H_1=2.0)) for t in (:DOM, :H_1, :NH₄)] ≈ [-3.0, 0.75, 2.25]

    selective_temperature = _controlled(; parameters=(temperature_q10=(P_1=2.0, P_2=1.0),))
    warm = _fields(; NO₃=1.0, T=30.0, P_1=1.0, P_2=1.0)
    p1 = selective_temperature(1, 1, 1, _GRID, Val(:P_1), (; time=0.0), warm, aux)
    p2 = selective_temperature(1, 1, 1, _GRID, Val(:P_2), (; time=0.0), warm, aux)
    @test p1 ≈ 2 * p2
end
