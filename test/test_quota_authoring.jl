using Test

using Agate.Components: Plankton, Pool
using Agate.Processes:
    FixedStoichiometry, Growth, Liebig, ModelDefinition, Monod, NormalizedDroop,
    NutrientResponse, NutrientLimitation, QuotaResponse


@testset "Quota structural errors fail during canonicalization" begin
    components = quota_components()
    light = quota_processes().growth.factors.light
    valid_nitrogen = quota_response(
        :nitrogen, :minimum_nitrogen_quota, :maximum_nitrogen_quota
    )
    valid_phosphorus = quota_response(
        :phosphorus, :minimum_phosphorus_quota, :maximum_phosphorus_quota
    )
    growth_with(
        responses; stoichiometry=nothing, additional_resources=NamedTuple()
    ) = Growth(;
        plankton=:P,
        reference_resource=:DIC,
        additional_resources=additional_resources,
        bindings=(maximum_rate=:maximum_growth_rate,),
        stoichiometry=stoichiometry,
        factors=(light=light, nutrients=NutrientLimitation(Liebig(); responses=responses)),
    )
    definition_with(processes; components=components) = ModelDefinition(;
        components=components, processes=processes, parameters=quota_parameters()
    )

    mixed = growth_with((
        nitrogen=NutrientResponse(
            Monod(); resource=:DIN, bindings=(half_saturation=:nitrogen_half_saturation,)
        ),
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
    conflicting_growth = growth_with(
        (nitrogen=valid_nitrogen, phosphorus=valid_phosphorus);
        additional_resources=(nitrogen=:DIN,),
        stoichiometry=FixedStoichiometry(; reference_element=:carbon),
    )
    reference_quota = growth_with((
        carbon=QuotaResponse(
            NormalizedDroop();
            variable_state=:carbon,
            bindings=(
                minimum_quota=:minimum_nitrogen_quota,
                maximum_quota=:maximum_nitrogen_quota,
            ),
        ),
    ))

    nitrogen_bindings = quota_processes().nitrogen_uptake.bindings
    wrong_element_uptake = quota_uptake(:nitrogen, :PO4, nitrogen_bindings)
    incomplete_uptake = quota_uptake(
        :nitrogen, :DIN, (maximum_rate=:maximum_nitrogen_uptake,)
    )

    mixed_parameters = (
        maximum_growth_rate=Agate.Parameters.Parameter(0.5),
        photosynthetic_slope=Agate.Parameters.Parameter(0.05),
        nitrogen_half_saturation=Agate.Parameters.Parameter(0.2),
        minimum_phosphorus_quota=Agate.Parameters.Parameter(0.005),
        maximum_phosphorus_quota=Agate.Parameters.Parameter(0.02),
    )
    @test_nowarn Agate.Processes.canonicalize_model(ModelDefinition(;
        components, processes=(growth=mixed,), parameters=mixed_parameters
    ))

    cases = (
        (definition_with((growth=unknown_state,)), ("no state :missing",)),
        (
            definition_with((growth=conflicting_growth,)),
            ("additional resource :nitrogen conflicts", "explicit prognostic state :nitrogen"),
        ),
        (
            definition_with((growth=reference_quota,)),
            ("quota response variable state :carbon", "must differ from reference state :carbon"),
        ),
        (definition_with((uptake=incomplete_uptake,)), ("requires explicit bindings",)),
        (
            definition_with((uptake=wrong_element_uptake,)),
            ("element :phosphorus", "expected :nitrogen"),
        ),
    )

    for (definition, fragments) in cases
        message = canonicalization_error_message(definition)
        @test all(fragment -> occursin(fragment, message), fragments)
    end

end
