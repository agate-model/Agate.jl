using Test
using Oceananigans.Architectures: CPU
using Oceananigans.Grids: RectilinearGrid
using Oceananigans.Fields: interior
using Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    DissolvedParticulate, LOBSTER, NutrientsPlanktonDetritus, nutrient_uptake,
    solid_waste, dissolved_waste, inorganic_waste
using OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients, NitrateAmmonia, N, PO₄

# A directly authored model exercises Agate's NPD adapter independently of providers.
function npd_test_definition()
    C = Agate.Components
    P = Agate.Processes
    components = (
        NO₃=C.Pool(:nitrogen), NH₄=C.Pool(:nitrogen), DOM=C.Pool(:nitrogen),
        solid_waste=C.Pool(:nitrogen), dissolved_waste=C.Pool(:nitrogen),
        inorganic_waste=C.Pool(:nitrogen),
        P=C.Plankton(;
            states=(nitrogen=:nitrogen,), reference_state=:nitrogen,
            size_structure=[1.0],
        ),
    )
    processes = (
        growth=P.Growth(;
            plankton=:P, reference_resource=:NO₃,
            bindings=(maximum_rate=:maximum_growth_rate,),
            factors=(
                light=P.Light(P.Smith(); driver=:PAR),
                temperature=P.Temperature(P.Q10(:plankton); driver=:T,
                    bindings=(q10=:temperature_q10, reference_temperature=:reference_temperature)),
            ),
        ),
        mortality=P.Mortality(P.LinearMortality();
            plankton=:P, bindings=(rate=:mortality_rate,),
            products=P.Products((
                solid=:solid_waste, dissolved=:dissolved_waste, inorganic=:inorganic_waste,
            ); fractions=(dissolved=:dissolved_fraction, inorganic=:inorganic_fraction)),
        ),
        remineralization=P.Remineralization(P.LinearRemineralization();
            sources=:DOM, destination=:NH₄, bindings=(rate=:remineralization_rate,),
        ),
    )
    parameters = (
        maximum_growth_rate=Agate.Parameters.Parameter(1.0),
        alpha=Agate.Parameters.Parameter(0.1),
        temperature_q10=Agate.Parameters.Parameter(2.0),
        reference_temperature=Agate.Parameters.Parameter(20.0),
        mortality_rate=Agate.Parameters.Parameter(0.2),
        dissolved_fraction=Agate.Parameters.Parameter(0.25),
        inorganic_fraction=Agate.Parameters.Parameter(0.125),
        remineralization_rate=Agate.Parameters.Parameter(0.4),
    )
    return Agate.ModelDefinition(; components, processes, parameters)
end

npd_test_cell(x) = fill(x, 1, 1, 1)
npd_test_fields(; NO₃=1.0, NH₄=0.1, DOM=0.5, P_1=1.0, T=20.0) = (
    NO₃=npd_test_cell(NO₃), NH₄=npd_test_cell(NH₄),
    DOM=npd_test_cell(DOM), sPOM=npd_test_cell(0.0), bPOM=npd_test_cell(0.0),
    P_1=npd_test_cell(P_1), T=npd_test_cell(T), PO₄=npd_test_cell(1.0),
)

# The external provider contract is exercised without depending on a model package.
struct TestNPDFamily <: Agate.ModelFamilies.AbstractModelFamily end
Agate.Construction.family_id(::TestNPDFamily) = :TestNPD
Agate.ModelFamilies.default_components(::TestNPDFamily) = npd_test_definition().components
Agate.ModelFamilies.default_processes(::TestNPDFamily) = npd_test_definition().processes
Agate.Parameters.parameter_definitions(::TestNPDFamily) = npd_test_definition().parameters
Agate.ModelFamilies.plankton_roles(::TestNPDFamily) = (phytoplankton=:P,)
Agate.ModelFamilies.setting_definitions(::TestNPDFamily) = (
    chlorophyll_ratio=Agate.ModelFamilies.ModelSetting(1.31; domain=:nonnegative),
    carbon_ratio=Agate.ModelFamilies.ModelSetting(6.56; domain=:positive),
)
Agate.Integrations.npd_configuration(::TestNPDFamily, _) = (
    nutrient_tracers=(:NO₃, :NH₄), consumed_detritus=(:DOM,), dependencies=(:T,),
)

@testset "NPD model-provider interface" begin
    C, NPD = Agate.Construction, Agate.Integrations
    family = TestNPDFamily()
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    pfts = C.plankton_realization(family, (phytoplankton=(P=[1.0],),))
    @test keys(pfts) == (:P,)
    @test Agate.ModelFamilies.plankton_roles(family).phytoplankton == :P
    plankton, recipe = NPD.construct_npd_plankton_plus_recipe(family;
        plankton_pfts=pfts, setting_overrides=(chlorophyll_ratio=1.5,), grid,
        scalar_type=Float32,
    )
    @test C.decode_recipe(C.encode_recipe(recipe)) == recipe
    replayed = NPD.construct_npd_plankton(family, C.decode_recipe(C.encode_recipe(recipe));
        grid, scalar_type=Float32)
    direct = NPD.construct_npd_plankton(family;
        plankton_pfts=pfts, setting_overrides=(chlorophyll_ratio=1.5,), grid,
        scalar_type=Float32)
    @test replayed.runtime.parameters == plankton.runtime.parameters == direct.runtime.parameters
    @test replayed.settings == plankton.settings == direct.settings
    @test plankton.settings === plankton.runtime.metadata.model_settings
    @test plankton.settings.chlorophyll_ratio === 1.5f0
    @test plankton.settings.carbon_ratio === Float32(6.56)
    @test Agate.Introspection.model_settings(plankton) == plankton.settings
    @test Agate.Integrations._exchange_tracers(plankton) ==
          (:solid_waste, :dissolved_waste, :inorganic_waste)
    @test required_biogeochemical_auxiliary_fields(plankton) == (:PAR,)
    @test :T in required_biogeochemical_tracers(LOBSTER(grid; plankton).underlying_biogeochemistry)

    config = (
        owned_components=(:P,), phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃, :NH₄), consumed_detritus=(:DOM,), dependencies=(:T,),
    )
    @test NPD.NPDPlankton(plankton.runtime; config...,
        settings=(carbon_ratio=6.56, chlorophyll_ratio=1.5)).settings === plankton.settings
    @test_throws "must match runtime.metadata.model_settings" NPD.NPDPlankton(
        plankton.runtime; config..., settings=(carbon_ratio=7.0, chlorophyll_ratio=1.5))
    @test_throws "Available settings" NPD.construct_npd_plankton(family;
        plankton_pfts=pfts, setting_overrides=(unknown=1.0,))
    @test_throws ArgumentError NPD.construct_npd_plankton(family;
        plankton_pfts=pfts, setting_overrides=(carbon_ratio=0.0,))
    @test_throws "grid is required" NPD.construct_npd_plankton(family;
        plankton_pfts=pfts, sinking_tracers=(P_1=0.1,))
    wrong = C.ModelRecipe(:Other, recipe.plankton_pfts, recipe.parameter_overrides,
        recipe.setting_overrides, recipe.sinking_tracers, recipe.open_bottom)
    @test_throws r"Recipe family" NPD.construct_npd_plankton(family, wrong)
end

@testset "Generic OceanBioME NPD coupling" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    runtime = Agate.Construction.construct(npd_test_definition(); grid)
    settings = (carbon_ratio=6.56, chlorophyll_ratio=1.31)
    base = (
        owned_components=(:P,), phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃, :NH₄), consumed_detritus=(:DOM,),
        dependencies=(:T,),
    )
    NPDPlankton = Agate.Integrations.NPDPlankton
    @test_throws "requires explicit settings" NPDPlankton(runtime; base...)
    plankton = NPDPlankton(runtime; base..., settings)
    @test required_biogeochemical_tracers(plankton) == (:P_1,)
    @test required_biogeochemical_auxiliary_fields(runtime) == (:PAR, :T)
    @test required_biogeochemical_auxiliary_fields(plankton) == (:PAR,)
    @test Agate.Introspection.model_settings(plankton) == settings
    @test Agate.Integrations._exchange_tracers(plankton) ==
          (:solid_waste, :dissolved_waste, :inorganic_waste)

    reordered = NPDPlankton(runtime; base..., settings, exchange_tracers=(
        inorganic=:inorganic_waste, solid=:solid_waste, dissolved=:dissolved_waste,
    ))
    @test typeof(reordered) === typeof(plankton)
    @test Agate.Integrations._exchange_tracers(reordered) ==
          Agate.Integrations._exchange_tracers(plankton)
    for invalid in (
        (; nutrient_tracers=(:NO₃,)),
        (; consumed_detritus=()),
        (; exchange_tracers=(solid=nothing, dissolved=:dissolved_waste,
                             inorganic=:inorganic_waste)),
    )
        @test_throws "no OceanBioME route" NPDPlankton(
            runtime; merge(base, invalid)..., settings
        )
    end
    for (invalid, message) in (
        ((; owned_components=(:P, :P)), "duplicate names"),
        ((; nutrient_tracers=(:NO₃, :NO₃)), "duplicate names"),
        ((; nutrient_tracers=(:P_1,)), "must not overlap"),
        ((; consumed_detritus=(:missing,)), "not an Agate runtime tracer"),
        ((; consumed_detritus=(:DOM,), nutrient_tracers=(:DOM,)), "must not overlap"),
        ((; exchange_tracers=(solid=:P_1, dissolved=nothing, inorganic=nothing)),
         "must not overlap"),
        ((; exchange_tracers=(solid=:solid_waste, dissolved=:solid_waste,
                              inorganic=nothing)), "must not overlap"),
        ((; dependencies=(:T, :T)), "duplicate names"),
        ((; owned_components=(:P,), phytoplankton_components=(:NO₃,)), "must belong"),
        ((; exchange_tracers=(solid=:solid_watse, dissolved=nothing,
                              inorganic=nothing)), "unknown runtime tracer"),
    )
        @test_throws message NPDPlankton(runtime; merge(base, invalid)..., settings)
    end
    for exchange_tracers in (
        (solid=:solid_waste, dissolved=:dissolved_waste),
        (solid=:solid_waste, dissolved=:dissolved_waste,
         inorganic=:inorganic_waste, extra=nothing),
    )
        @test_throws "must define exactly" NPDPlankton(runtime; base..., settings,
            exchange_tracers)
    end

    bgc = LOBSTER(grid; plankton).underlying_biogeochemistry
    @test all(name -> name in required_biogeochemical_tracers(bgc),
        (:NO₃, :NH₄, :DOM, :sPOM, :bPOM, :T, :P_1))
    make_npd(; nutrients=bgc.nutrients, detritus=bgc.detritus, plankton=bgc.plankton) =
        NutrientsPlanktonDetritus{Float64}(
            nutrients, plankton, detritus, bgc.inorganic_carbon, bgc.oxygen
        )
    @test required_biogeochemical_tracers(make_npd(; plankton=reordered)) ==
          required_biogeochemical_tracers(bgc)
    @test_throws ArgumentError required_biogeochemical_tracers(
        make_npd(; nutrients=Nutrients(; nitrogen=N))
    )
    @test_throws "consumed_detritus" required_biogeochemical_tracers(
        make_npd(; detritus=nothing)
    )
    @test_throws "also belong to native" required_biogeochemical_tracers(
        make_npd(; plankton=NPDPlankton(runtime;
            merge(base, (; owned_components=(:P, :DOM), consumed_detritus=()))..., settings))
    )

    fields = npd_test_fields()
    aux = (PAR=npd_test_cell(100.0),)
    uptake(nutrient, npd=bgc) = nutrient_uptake(
        1, 1, 1, grid, Val(nutrient), npd.plankton, npd, fields, aux
    )
    @test uptake(:NO₃) > 0
    @test uptake(:DOM) == 0
    @test solid_waste(1, 1, 1, grid, plankton, bgc, fields, aux) ≈ 0.125
    @test dissolved_waste(1, 1, 1, grid, plankton, bgc, fields, aux) ≈ 0.05
    @test inorganic_waste(1, 1, 1, grid, plankton, bgc, fields, aux) ≈ 0.025
    @test isfinite(bgc(1, 1, 1, grid, Val(:P_1), (; time=0.0), fields, aux))
    phosphate = make_npd(; nutrients=Nutrients(; nitrogen=NitrateAmmonia(), phosphate=PO₄))
    @test uptake(:PO₄, phosphate) ≈ (uptake(:NO₃) + uptake(:NH₄)) / 16
    @test uptake(:C) ≈ settings.carbon_ratio * (uptake(:NO₃) + uptake(:NH₄))

    # Model-independent sinking forwarding at the composed NPD boundary.
    column_grid = RectilinearGrid(CPU(); size=(1, 1, 4), extent=(1, 1, 40))
    sinking = Agate.Construction.construct(
        npd_test_definition(); grid=column_grid,
        sinking_tracers=(P_1=0.1,), open_bottom=true,
    )
    sinking_plankton = NPDPlankton(sinking; base..., settings)
    column = LOBSTER(column_grid; plankton=sinking_plankton).underlying_biogeochemistry
    drift = biogeochemical_drift_velocity(column, Val(:P_1)).w
    @test drift === biogeochemical_drift_velocity(sinking_plankton, Val(:P_1)).w
    @test any(v -> v < 0, Array(interior(drift)))
end

@testset "Generic NPD nutrient subsets and conservation" begin
    grid = RectilinearGrid(CPU(); size=(1, 1, 1), extent=(1, 1, 1))
    NPDPlankton = Agate.Integrations.NPDPlankton
    settings = (carbon_ratio=6.56, chlorophyll_ratio=1.31)

    # Native NH4 remains available but is not consumed by an NO3-only Agate model.
    definition = Agate.ModelDefinition(;
        components=(
            NO₃=Agate.Components.Pool(:nitrogen),
            P=Agate.Components.Plankton(;
                states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure=[1.0]),
        ),
        processes=(growth=Agate.Processes.Growth(;
            plankton=:P, reference_resource=:NO₃,
            bindings=(maximum_rate=:maximum_growth_rate,)),),
        parameters=(maximum_growth_rate=Agate.Parameters.Parameter(1.0),),
    )
    single = NPDPlankton(Agate.Construction.construct(definition; grid);
        owned_components=(:P,), phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃,), settings)
    npd = LOBSTER(grid; plankton=single).underlying_biogeochemistry
    fields = npd_test_fields(; NO₃=1.0, NH₄=1.0, P_1=1.0)
    uptake(nutrient) = nutrient_uptake(
        1, 1, 1, grid, Val(nutrient), single, npd, fields, (;))
    @test :NH₄ in required_biogeochemical_tracers(npd)
    @test uptake(:NO₃) ≈ 1.0
    @test uptake(:NH₄) == 0
    @test npd(1, 1, 1, grid, Val(:P_1), (; time=0.0), fields, (;)) ≈ 1.0

    # Unused exchange channels are absent in a model with no waste pools.
    bare = runtime_npz_model(; grid)
    wrapper = NPDPlankton(bare;
        owned_components=(:P, :Z), phytoplankton_components=(:P,),
        nutrient_tracers=(:N, :D), settings)
    @test Agate.Integrations._exchange_tracers(wrapper) == (nothing, nothing, nothing)

    # The composed OceanBioME boundary conserves nitrogen across all hooks.
    runtime = Agate.Construction.construct(npd_test_definition(); grid)
    coupled = NPDPlankton(runtime;
        owned_components=(:P,), phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃, :NH₄), consumed_detritus=(:DOM,),
        dependencies=(:T,), settings)
    detritus = DissolvedParticulate(grid;
        dissolved_remineralisation_rate=0.0,
        particulate_remineralisation_rate=(0.0, 0.0), sinking_speeds=(0.0, 0.0))
    closed = LOBSTER(grid; plankton=coupled, detritus,
        nutrients=Nutrients(; nitrogen=NitrateAmmonia(; nitrification_rate=0.0)),
    ).underlying_biogeochemistry
    state = npd_test_fields(; NO₃=2.0, NH₄=0.3, DOM=0.8, P_1=0.12)
    aux = (PAR=npd_test_cell(100.0),)
    tracers = (:NO₃, :NH₄, :DOM, :sPOM, :bPOM, :P_1)
    fluxes = map(tracers) do tracer
        closed(1, 1, 1, grid, Val(tracer), (; time=0.0), state, aux)
    end
    @test all(isfinite, fluxes)
    @test isapprox(sum(fluxes), 0; atol=1e-12)
end
