using Test

using Agate.Configuration: Population, Pool, population_state, realize_model_layout
using Agate.Parameters: Parameter
using Agate.Processes:
    FixedStoichiometry, FrankTNorm, Growth, Liebig, Light, ModelDefinition, Monod, NormalizedDroop,
    NutrientResponse, Nutrients, NutrientUptake, QuotaRegulatedMonod, QuotaResponse,
    Smith, build_parameter_plan, driver_identities, normalize_model, planned_parameter_slot


quota_components() = (
    DIC=Pool(:carbon),
    DIN=Pool(:nitrogen),
    PO4=Pool(:phosphorus),
    P=Population(;
        states=(carbon=:carbon, nitrogen=:nitrogen, phosphorus=:phosphorus),
        size_structure=[1.0, 2.0],
    ),
)

quota_response(state, minimum, maximum) = QuotaResponse(
    NormalizedDroop();
    target=population_state(:P, state),
    reference=population_state(:P, :carbon),
    bindings=(minimum_quota=minimum, maximum_quota=maximum),
)

quota_uptake(state, resource, bindings) = NutrientUptake(
    QuotaRegulatedMonod();
    population=:P,
    target_state=state,
    reference_state=:carbon,
    resource=resource,
    bindings=bindings,
)

function quota_processes()
    responses = (
        nitrogen=quota_response(
            :nitrogen, :minimum_nitrogen_quota, :maximum_nitrogen_quota
        ),
        phosphorus=quota_response(
            :phosphorus, :minimum_phosphorus_quota, :maximum_phosphorus_quota
        ),
    )
    growth = Growth(;
        populations=:P,
        state=:carbon,
        source=:DIC,
        factors=(
            light=Light(
                Smith(); driver=:PAR,
                bindings=(maximum_rate=:maximum_growth_rate, alpha=:photosynthetic_slope),
            ),
            nutrients=Nutrients(Liebig(); responses=responses),
        ),
    )
    nitrogen_uptake = quota_uptake(:nitrogen, :DIN, (
        maximum_rate=:maximum_nitrogen_uptake,
        K=:nitrogen_half_saturation,
        minimum_quota=:minimum_nitrogen_quota,
        maximum_quota=:maximum_nitrogen_quota,
        hill=:nitrogen_uptake_hill,
    ))
    phosphorus_uptake = quota_uptake(:phosphorus, :PO4, (
        maximum_rate=:maximum_phosphorus_uptake,
        K=:phosphorus_half_saturation,
        minimum_quota=:minimum_phosphorus_quota,
        maximum_quota=:maximum_phosphorus_quota,
        hill=:phosphorus_uptake_hill,
    ))
    return (; growth, nitrogen_uptake, phosphorus_uptake)
end

quota_parameters() = (
    maximum_growth_rate=Parameter(0.5; axes=:plankton),
    photosynthetic_slope=Parameter(0.05; axes=:plankton),
    minimum_nitrogen_quota=Parameter(0.05; axes=:plankton),
    maximum_nitrogen_quota=Parameter(0.2; axes=:plankton),
    minimum_phosphorus_quota=Parameter(0.005; axes=:plankton),
    maximum_phosphorus_quota=Parameter(0.02; axes=:plankton),
    maximum_nitrogen_uptake=Parameter(0.1; axes=:plankton),
    nitrogen_half_saturation=Parameter(0.2; axes=:plankton),
    nitrogen_uptake_hill=Parameter(2.0; axes=:plankton),
    maximum_phosphorus_uptake=Parameter(0.01; axes=:plankton),
    phosphorus_half_saturation=Parameter(0.02; axes=:plankton),
    phosphorus_uptake_hill=Parameter(2.0; axes=:plankton),
)

quota_definition() = ModelDefinition(;
    components=quota_components(), processes=quota_processes(), parameters=quota_parameters()
)

function quota_normalization_error(definition)
    error = try
        normalize_model(definition)
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    return error isa Exception ? sprint(showerror, error) : ""
end

@testset "Droop quota authoring and normalization" begin
    definition = quota_definition()
    normalized = normalize_model(definition)
    growth = normalized.processes.growth

    @test growth.facts.routing == (mode=:quota, factor=:nutrients, source=:DIC)
    @test growth.facts.population_states == (population_state(:P, :carbon),)
    @test normalized.processes.nitrogen_uptake.facts == (
        target=population_state(:P, :nitrogen),
        reference=population_state(:P, :carbon),
        resource=:DIN,
    )
    @test driver_identities(normalized) == (:PAR,)
    @test Nutrients(
        FrankTNorm(); responses=(nitrogen=quota_response(
            :nitrogen, :minimum_nitrogen_quota, :maximum_nitrogen_quota
        ),)
    ) isa Nutrients

    explicit_state = normalize_model(ModelDefinition(;
        components=(
            P=Population(; states=(carbon=:carbon, nitrogen=:nitrogen)),
            DIC=Pool(:carbon),
        ),
        processes=(growth=Growth(;
            populations=:P,
            state=:carbon,
            factors=(
                light=Light(Smith(); driver=:PAR),
                nutrients=NutrientResponse(Monod(); resource=:DIC),
            ),
        ),),
    ))
    @test explicit_state.processes.growth.facts.population_states ==
        (population_state(:P, :carbon),)

    layout = realize_model_layout(
        definition.components;
        scalar_type=Float64,
        auxiliary_fields=driver_identities(normalized),
    )
    plan = build_parameter_plan(normalized, layout)
    nitrogen_quota = plan.parameters.minimum_nitrogen_quota
    @test (nitrogen_quota.storage_shape, nitrogen_quota.storage_labels) == (
        (2,), ((:P_1, :P_2),)
    )
    @test nitrogen_quota.applicable_indices == ((1, 2),)

    shared = Tuple(
        binding for binding in normalized.parameter_bindings
        if binding.parameter === :minimum_nitrogen_quota
    )
    @test length(shared) == 2
    @test all(planned_parameter_slot(plan, binding) == ((1, 2),) for binding in shared)
end

@testset "Quota structural errors fail during normalization" begin
    components = quota_components()
    light = quota_processes().growth.factors.light
    valid_nitrogen = quota_response(
        :nitrogen, :minimum_nitrogen_quota, :maximum_nitrogen_quota
    )
    valid_phosphorus = quota_response(
        :phosphorus, :minimum_phosphorus_quota, :maximum_phosphorus_quota
    )
    growth_with(responses; state=:carbon, stoichiometry=nothing) = Growth(;
        populations=:P,
        state=state,
        source=:DIC,
        stoichiometry=stoichiometry,
        factors=(light=light, nutrients=Nutrients(Liebig(); responses=responses)),
    )
    definition_with(processes; components=components) = ModelDefinition(;
        components=components, processes=processes, parameters=quota_parameters()
    )

    mixed = growth_with((
        nitrogen=NutrientResponse(Monod(); resource=:DIN),
        phosphorus=valid_phosphorus,
    ))
    unknown_state = growth_with((
        nitrogen=QuotaResponse(
            NormalizedDroop();
            target=population_state(:P, :missing),
            reference=population_state(:P, :carbon),
            bindings=(
                minimum_quota=:minimum_nitrogen_quota,
                maximum_quota=:maximum_nitrogen_quota,
            ),
        ),
        phosphorus=valid_phosphorus,
    ))
    mismatched_population = growth_with((
        nitrogen=QuotaResponse(
            NormalizedDroop();
            target=population_state(:P, :nitrogen),
            reference=population_state(:Q, :carbon),
            bindings=(
                minimum_quota=:minimum_nitrogen_quota,
                maximum_quota=:maximum_nitrogen_quota,
            ),
        ),
        phosphorus=valid_phosphorus,
    ))
    missing_growth_state = growth_with(
        (nitrogen=valid_nitrogen, phosphorus=valid_phosphorus); state=nothing
    )
    fixed_quota_growth = growth_with(
        (nitrogen=valid_nitrogen, phosphorus=valid_phosphorus);
        stoichiometry=FixedStoichiometry(; reference=:carbon),
    )

    nitrogen_bindings = quota_processes().nitrogen_uptake.bindings
    wrong_currency_uptake = quota_uptake(:nitrogen, :PO4, nitrogen_bindings)
    structured_uptake = quota_uptake(:nitrogen, :DIN, nitrogen_bindings)
    incomplete_uptake = quota_uptake(
        :nitrogen, :DIN, (maximum_rate=:maximum_nitrogen_uptake,)
    )

    cases = (
        (definition_with((growth=mixed,)), ("all NutrientResponse", "all QuotaResponse")),
        (definition_with((growth=unknown_state,)), ("no state :missing",)),
        (
            definition_with((growth=mismatched_population,); components=merge(components, (
                Q=Population(; states=(carbon=:carbon,), size_structure=[1.0, 2.0]),
            ))),
            ("growth population :P",),
        ),
        (definition_with((growth=missing_growth_state,)), ("explicit state selection",)),
        (definition_with((growth=fixed_quota_growth,)), ("independent NutrientUptake",)),
        (definition_with((uptake=incomplete_uptake,)), ("requires explicit bindings",)),
        (
            definition_with((uptake=wrong_currency_uptake,)),
            ("currency :phosphorus", "expected :nitrogen"),
        ),
        (
            definition_with(
                (uptake=structured_uptake,);
                components=merge(components, (DIN=Pool(:nitrogen; size_structure=[1.0, 2.0]),)),
            ),
            ("nutrient uptake resource", "scalar Pool"),
        ),
    )

    for (definition, fragments) in cases
        message = quota_normalization_error(definition)
        @test all(fragment -> occursin(fragment, message), fragments)
    end
end
