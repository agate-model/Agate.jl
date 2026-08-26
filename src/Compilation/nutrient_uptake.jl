function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ModelLayout,
    plan::ParameterPlan,
) where {P<:NutrientUptake}
    process = named.process
    facts = named.facts
    target_tracers, _ = _realize_population_state(facts.target, layout)
    reference_tracers, _ = _realize_population_state(facts.reference, layout)
    resource = _scalar_component_target(layout, facts.resource)
    slots = parameter_slot_bindings(definition, named, (), process)
    fluxes = Any[]

    for population_axis in eachindex(target_tracers)
        axis_positions = (population=_axis_position(population_axis),)
        operands = (
            input_operand(layout, resource),
            input_operand(layout, target_tracers[population_axis]),
            input_operand(layout, reference_tracers[population_axis]),
            parameter_operand(slots.maximum_rate, plan, axis_positions),
            parameter_operand(slots.K, plan, axis_positions),
            parameter_operand(slots.minimum_quota, plan, axis_positions),
            parameter_operand(slots.maximum_quota, plan, axis_positions),
            parameter_operand(slots.hill, plan, axis_positions),
        )
        rate = RateElement(formulation(process), operands)
        push!(
            fluxes,
            FluxSpec(resource, rate, Weight{-1}()),
            FluxSpec(target_tracers[population_axis], rate, Weight{1}()),
        )
    end
    return Tuple(fluxes)
end
