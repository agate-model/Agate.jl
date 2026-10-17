using Test
using Oceananigans.Architectures: CPU
using Oceananigans.Grids: RectilinearGrid
using Oceananigans.Fields: ConstantField, ZeroField, interior
using Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers

using OceanBioME: chlorophyll, PrescribedPhotosyntheticallyActiveRadiation
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    DissolvedParticulate, LOBSTER, NutrientsPlanktonDetritus, nutrient_uptake,
    solid_waste, dissolved_waste, inorganic_waste
using OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients, NitrateAmmonia, N, PO₄

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
        grid=_GRID, arch=CPU(), scalar_type=Float32, size_structure, parameters, settings,
        sinking_tracers=(nano_1=0.1,), open_bottom=false,
    )
    encoded = Agate.Construction.encode_recipe(recipe)
    decoded = Agate.Construction.decode_recipe(encoded)
    replayed = FrankenLOBSTER.construct(decoded; grid=_GRID, arch=CPU(), scalar_type=Float32)
    direct32 = FrankenLOBSTER.construct(; arch=CPU(), scalar_type=Float32)
    bgc = LOBSTER(_GRID; plankton)

    # The NPD composite must forward plankton sinking to Agate, while native
    # OceanBioME detritus drift retains its own dispatch.
    npd = bgc.underlying_biogeochemistry
    sinking = biogeochemical_drift_velocity(npd, Val(:nano_1)).w
    @test sinking === biogeochemical_drift_velocity(plankton, Val(:nano_1)).w
    @test sinking === biogeochemical_drift_velocity(bgc, Val(:nano_1)).w
    @test biogeochemical_drift_velocity(npd, Val(:pico_1)).w == ZeroField()
    @test biogeochemical_drift_velocity(npd, Val(:sPOM)).w ===
          biogeochemical_drift_velocity(bgc, Val(:sPOM)).w

    # The single-layer grid above has a closed bottom, which can suppress
    # sinking at every face. Test nonzero transport in a resolved open column.
    column_grid = RectilinearGrid(CPU(), Float32; size=(1, 1, 4), extent=(1f0, 1f0, 40f0))
    column_plankton = FrankenLOBSTER.construct(;
        grid=column_grid, sinking_tracers=(P_1=0.1f0,), open_bottom=true
    )
    column_npd = LOBSTER(column_grid; plankton=column_plankton).underlying_biogeochemistry
    column_drift = biogeochemical_drift_velocity(column_npd, Val(:P_1)).w
    @test any(v -> v < 0, Array(interior(column_drift)))

    @test required_biogeochemical_tracers(plankton) ==
          (:nano_1, :pico_1, :meso_1, :micro_1, :heterotroph_1, :heterotroph_2)
    @test size(plankton.runtime.parameters.bacterial_dom_half_saturation) == (2, 1)
    @test :T ∉ required_biogeochemical_tracers(plankton.runtime)
    @test required_biogeochemical_auxiliary_fields(plankton.runtime) == (:PAR, :T)
    @test required_biogeochemical_auxiliary_fields(plankton) == (:PAR,)
    @test all(t -> t in required_biogeochemical_tracers(bgc), (:NO₃, :NH₄, :DOM, :sPOM, :bPOM, :T))
    @test chlorophyll(plankton, (tracers=(nano_1=_cell(2.0), pico_1=_cell(1.0)),))[1, 1, 1] ≈ 4.5
    @test recipe.setting_overrides == settings
    @test encoded["family"] == "FrankenLOBSTER"
    @test encoded["realization"]["setting_overrides"]["chlorophyll_ratio"] == 1.5
    mktemp() do path, io
        close(io)
        Agate.Construction.export_recipe(path, recipe)
        @test Agate.Construction.import_recipe(path) == recipe
    end

    I = Agate.Introspection
    for inspect in (
        I.model_settings, I.pfts, I.parameter_names, I.plankton_tracers,
        I.plankton_diameters, I.model_summary,
    )
        @test inspect(plankton) == inspect(plankton.runtime)
    end
    @test plankton.settings == I.model_settings(plankton)
    @test plankton.settings === plankton.runtime.metadata.model_settings
    @test I.parameter_domains(plankton, :assimilation_matrix) ==
          I.parameter_domains(plankton.runtime, :assimilation_matrix)
    @test I.interaction_matrix(plankton, :assimilation_matrix) ==
          I.interaction_matrix(plankton.runtime, :assimilation_matrix)

    @test (replayed.runtime.parameters, replayed.settings) == (plankton.runtime.parameters, plankton.settings)
    @test eltype(plankton.runtime.parameters.maximum_growth_rate) === Float32
    @test eltype(direct32.runtime.parameters.maximum_growth_rate) === Float32
    @test_throws ArgumentError Agate.Integrations.NPDPlankton(
        plankton.runtime;
        owned_components=(:P, :Z, :H),
        phytoplankton_components=(:P,),
        exchange_tracers=(solid=:solid_watse, dissolved=nothing, inorganic=nothing),
        settings=plankton.settings,
    )
    # The direct adapter accepts arbitrary nutrient names, but every nonzero
    # runtime equation must be owned or connected to exactly one external hook.
    base = (
        owned_components=(:P, :Z, :H),
        phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃, :NH₄),
        consumed_detritus=(:DOM,),
    )
    default_settings = Agate.Integrations.NPDPlankton(plankton.runtime; base...)
    @test default_settings.settings === plankton.runtime.metadata.model_settings
    @test Agate.Integrations.NPDPlankton(
        plankton.runtime; base..., settings=plankton.settings
    ).settings === plankton.settings
    @test I.model_settings(default_settings) == default_settings.settings
    # The explicit low-level keyword remains valid, including Float64 inputs to
    # a Float32 runtime, while incompatible ratios are rejected.
    @test Agate.Integrations.NPDPlankton(
        plankton.runtime;
        merge(base, (; settings=(carbon_ratio=6.56, chlorophyll_ratio=1.5)))...
    ).settings === plankton.settings
    for ratios in (
        (carbon_ratio=7.0, chlorophyll_ratio=1.5),
        (carbon_ratio=6.56, chlorophyll_ratio=1.6),
    )
        @test_throws "must match runtime.metadata.model_settings" Agate.Integrations.NPDPlankton(
            plankton.runtime; merge(base, (; settings=ratios))...
        )
    end
    canonical = Agate.Integrations.NPDPlankton(plankton.runtime; base...)
    reordered = Agate.Integrations.NPDPlankton(
        plankton.runtime;
        merge(base, (; exchange_tracers=(
            inorganic=:inorganic_waste, solid=:solid_waste, dissolved=:dissolved_waste,
        )))...
    )
    @test typeof(reordered) === typeof(canonical)
    @test Agate.Integrations._exchange_tracers(reordered) ==
          Agate.Integrations._exchange_tracers(canonical)
    for channels in (
        (solid=:solid_waste, dissolved=:dissolved_waste),
        (solid=:solid_waste, dissolved=:dissolved_waste, inorganic=:inorganic_waste, extra=nothing),
    )
        @test_throws "must define exactly" Agate.Integrations.NPDPlankton(
            plankton.runtime; merge(base, (; exchange_tracers=channels))...
        )
    end
    @test Agate.Integrations.NPDPlankton(
        plankton.runtime;
        merge(base, (; nutrient_tracers=(:NO₃, :NH₄, :DOM), consumed_detritus=()))...
    ) isa Agate.Integrations.NPDPlankton
    for missing_route in (
        (; nutrient_tracers=(:NO₃,)),
        (; consumed_detritus=()),
        (; exchange_tracers=(solid=nothing, dissolved=:dissolved_waste,
                              inorganic=:inorganic_waste)),
    )
        @test_throws "no OceanBioME route" Agate.Integrations.NPDPlankton(
            plankton.runtime; merge(base, missing_route)...
        )
    end
    for (invalid, message) in (
        ((; owned_components=(:P, :P, :Z, :H)), "duplicate names"),
        ((; owned_components=(:P, :Z), phytoplankton_components=(:H,)), "must belong"),
        ((; nutrient_tracers=(:nano_1,)), "must not overlap"),
        ((; nutrient_tracers=(:NO₃, :NO₃)), "duplicate names"),
        ((; consumed_detritus=(:missing,)), "not an Agate runtime tracer"),
        ((; consumed_detritus=(:DOM,), nutrient_tracers=(:DOM,)), "must not overlap"),
        ((; exchange_tracers=(solid=:nano_1, dissolved=nothing, inorganic=nothing)), "must not overlap"),
        ((; exchange_tracers=(solid=:solid_waste, dissolved=:solid_waste,
                               inorganic=nothing)), "must not overlap"),
        ((; dependencies=(:T, :T)), "duplicate names"),
    )
        @test_throws message Agate.Integrations.NPDPlankton(
            plankton.runtime; merge(base, invalid)...
        )
    end
    @test_throws "Available settings: chlorophyll_ratio, carbon_ratio" FrankenLOBSTER.construct(
        settings=(unknown=1.0,)
    )
    @test_throws ArgumentError FrankenLOBSTER.construct(sinking_tracers=(P_1=0.1,))

    # Generic Agate runtimes can still supply conversion ratios explicitly.
    # Their absent waste channels default to `nothing` rather than unknown names.
    nipizd = Agate.Models.NiPiZD.construct(; grid=_GRID)
    generic = (
        owned_components=(:P, :Z), phytoplankton_components=(:P,),
        nutrient_tracers=(:N, :D),
    )
    @test isempty(nipizd.metadata.model_settings)
    @test_throws "requires explicit settings" Agate.Integrations.NPDPlankton(
        nipizd; generic...
    )
    generic_settings = (carbon_ratio=6.56, chlorophyll_ratio=1.31)
    wrapper = Agate.Integrations.NPDPlankton(
        nipizd; generic..., settings=generic_settings
    )
    @test wrapper.settings == generic_settings
    @test I.model_settings(wrapper) == generic_settings
    @test Agate.Integrations._exchange_tracers(wrapper) == (nothing, nothing, nothing)
end

@testset "FrankenLOBSTER LOBSTER exchange" begin
    bgc = _controlled()
    # Validate configured exchange routes against the native NPD component
    # tracers before Oceananigans sets up its prognostic fields.
    make_npd(; nutrients=bgc.nutrients, detritus=bgc.detritus, plankton=bgc.plankton) =
        NutrientsPlanktonDetritus{Float64}(
            nutrients, plankton, detritus, bgc.inorganic_carbon, bgc.oxygen
        )
    @test required_biogeochemical_tracers(make_npd()) == required_biogeochemical_tracers(bgc)
    reordered = Agate.Integrations.NPDPlankton(
        bgc.plankton.runtime;
        owned_components=(:P, :Z, :H), phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃, :NH₄), consumed_detritus=(:DOM,),
        dependencies=(:T,),
        exchange_tracers=(inorganic=:inorganic_waste,
                          dissolved=:dissolved_waste, solid=:solid_waste),
    )
    @test required_biogeochemical_tracers(make_npd(; plankton=reordered)) ==
          required_biogeochemical_tracers(make_npd())
    @test_throws "nutrient_tracers" required_biogeochemical_tracers(
        make_npd(; nutrients=Nutrients(; nitrogen=N))
    )
    @test_throws "nutrient_tracers" required_biogeochemical_tracers(
        make_npd(; nutrients=Nutrients(; nitrogen=NitrateAmmonia(), phosphate=PO₄))
    )
    @test_throws "consumed_detritus" required_biogeochemical_tracers(
        make_npd(; detritus=nothing)
    )
    owns_dom = Agate.Integrations.NPDPlankton(
        bgc.plankton.runtime;
        owned_components=(:P, :Z, :H, :DOM),
        phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃, :NH₄),
        exchange_tracers=(solid=:solid_waste, dissolved=:dissolved_waste,
                          inorganic=:inorganic_waste),
        settings=bgc.plankton.settings,
    )
    @test_throws "also belong to native" required_biogeochemical_tracers(
        make_npd(; plankton=owns_dom)
    )
    aux = (PAR=_cell(1.0),)
    exchanged_state = _fields(; NO₃=1.0, DOM=0.3, P_1=2.0, Z_1=1.0)
    for hook in (solid_waste, dissolved_waste, inorganic_waste)
        @test hook(1, 1, 1, _GRID, reordered, make_npd(; plankton=reordered), exchanged_state, aux) ≈
              hook(1, 1, 1, _GRID, bgc.plankton, bgc, exchanged_state, aux)
    end
    tendency(tracer, fields) = bgc(1, 1, 1, _GRID, Val(tracer), (; time=0.0), fields, aux)
    uptake(tracer, fields) = nutrient_uptake(1, 1, 1, _GRID, Val(tracer), bgc.plankton, bgc, fields, aux)

    @test uptake(:DOM, _fields()) == 0
    # OceanBioME's elemental Val methods must remain unambiguous and retain
    # stoichiometric demand when the element lacks an explicit Agate equation.
    fallback_state = _fields(; NO₃=1.0, P_1=2.0)
    total_n_uptake = uptake(:NO₃, fallback_state) + uptake(:NH₄, fallback_state)
    for (nutrient, ratio) in (
        (:N, 1.0), (:PO₄, 1 / 16), (:Fe, 0.0032 / 16),
        (:Si, 0.0), (:C, bgc.plankton.settings.carbon_ratio),
    )
        @test nutrient_uptake(
            1, 1, 1, _GRID, Val(nutrient), bgc.plankton, bgc, fallback_state, aux
        ) ≈ ratio * total_n_uptake
    end

    state = _fields(; NO₃=1.0, NH₄=1.0, DOM=1.0, P_1=1.0, Z_1=1.0)
    @test all(
        tracer -> isfinite(tendency(tracer, state)),
        (:P_1, :P_2, :Z_1, :Z_2, :H_1, :NO₃, :NH₄, :DOM, :sPOM, :bPOM),
    )

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

@testset "FrankenLOBSTER closed nitrogen budget" begin
    bgc = _controlled(; parameters=(
        maximum_predation_rate=(Z_1=1.0, Z_2=0.5),
        phytoplankton_mortality_rate=(P_1=0.01, P_2=0.02),
        zooplankton_mortality_rate=(Z_1=0.01, Z_2=0.02),
        bacterioplankton_mortality_rate=(H_1=0.01,),
    ))
    state = _fields(;
        NO₃=2.0, NH₄=0.3, DOM=0.8, sPOM=0.2, bPOM=0.1,
        P_1=0.12, P_2=0.05, Z_1=0.06, Z_2=0.02, H_1=0.08,
    )
    nitrogen_tracers = (:NO₃, :NH₄, :DOM, :sPOM, :bPOM, :P_1, :P_2, :Z_1, :Z_2, :H_1)
    aux = (PAR=_cell(1.0),)
    fluxes = map(nitrogen_tracers) do tracer
        bgc(1, 1, 1, _GRID, Val(tracer), (; time=0.0), state, aux)
    end
    @test all(isfinite, fluxes)
    @test isapprox(sum(fluxes), 0.0; atol=1e-12)
end
