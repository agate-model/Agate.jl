using Test
using Oceananigans.Architectures: CPU
using Oceananigans.Grids: RectilinearGrid
using Oceananigans.Fields: interior
using Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    LOBSTER, NutrientsPlanktonDetritus, nutrient_uptake,
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
            products=P.Products(:solid_waste),
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
    )
        @test_throws message NPDPlankton(runtime; merge(base, invalid)..., settings)
    end
    @test_throws "must define exactly" NPDPlankton(runtime; base..., settings,
        exchange_tracers=(solid=:solid_waste, dissolved=:dissolved_waste))

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
    @test solid_waste(1, 1, 1, grid, plankton, bgc, fields, aux) ≈ 0.2
    @test dissolved_waste(1, 1, 1, grid, plankton, bgc, fields, aux) == 0
    @test inorganic_waste(1, 1, 1, grid, plankton, bgc, fields, aux) == 0
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
