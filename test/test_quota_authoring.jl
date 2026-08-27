using Test

using Agate.Configuration: Population, Pool, population_state
using Agate.Parameters: Parameter
using Agate.Processes:
    FixedStoichiometry, Growth, Liebig, Light, ModelDefinition, Monod, NormalizedDroop,
    NutrientResponse, Nutrients, NutrientUptake, QuotaRegulatedMonod, QuotaResponse, Smith


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
        bindings=(maximum_rate=:maximum_growth_rate,),
        state=:carbon,
        source=:DIC,
        factors=(
            light=Light(
                Smith(); driver=:PAR,
                bindings=(alpha=:photosynthetic_slope,),
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
    maximum_growth_rate=Parameter(0.5),
    photosynthetic_slope=Parameter(0.05),
    minimum_nitrogen_quota=Parameter(0.05),
    maximum_nitrogen_quota=Parameter(0.2),
    minimum_phosphorus_quota=Parameter(0.005),
    maximum_phosphorus_quota=Parameter(0.02),
    maximum_nitrogen_uptake=Parameter(0.1),
    nitrogen_half_saturation=Parameter(0.2),
    nitrogen_uptake_hill=Parameter(2.0),
    maximum_phosphorus_uptake=Parameter(0.01),
    phosphorus_half_saturation=Parameter(0.02),
    phosphorus_uptake_hill=Parameter(2.0),
)

quota_definition() = ModelDefinition(;
    components=quota_components(), processes=quota_processes(), parameters=quota_parameters()
)

@testset "Quota structural errors fail during canonicalization" begin
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
        bindings=(maximum_rate=:maximum_growth_rate,),
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
        message = canonicalization_error_message(definition)
        @test all(fragment -> occursin(fragment, message), fragments)
    end
end
