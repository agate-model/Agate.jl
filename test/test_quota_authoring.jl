using Test

using Agate.Configuration: Population, Pool
using Agate.Processes:
    FixedStoichiometry, Growth, Liebig, ModelDefinition, Monod, NormalizedDroop,
    NutrientResponse, Nutrients, QuotaResponse


@testset "Quota structural errors fail during canonicalization" begin
    components = quota_components()
    light = quota_processes().growth.factors.light
    valid_nitrogen = quota_response(
        :nitrogen, :minimum_nitrogen_quota, :maximum_nitrogen_quota
    )
    valid_phosphorus = quota_response(
        :phosphorus, :minimum_phosphorus_quota, :maximum_phosphorus_quota
    )
    growth_with(responses; stoichiometry=nothing) = Growth(;
        populations=:P,
        bindings=(maximum_rate=:maximum_growth_rate,),
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
            variable_state=:missing,
            bindings=(
                minimum_quota=:minimum_nitrogen_quota,
                maximum_quota=:maximum_nitrogen_quota,
            ),
        ),
        phosphorus=valid_phosphorus,
    ))
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

    canonical = Agate.Processes.canonicalize_model(quota_definition())
    @test only(canonical.processes.growth.facts.population_states).state === :carbon
    @test canonical.processes.nitrogen_uptake.facts.reference.state === :carbon
    @test canonical.processes.growth.process.factors.nutrients.responses.nitrogen.variable_state === :nitrogen
end
