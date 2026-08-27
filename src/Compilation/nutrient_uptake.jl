function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:NutrientUptake}
    process = named.process
    facts = named.facts
    layout = context.layout
    plan = context.plan
    target_tracers = state_tracers(layout, facts.target)
    reference_tracers = state_tracers(layout, facts.reference)
    resource = _scalar_component_target(layout, facts.resource)
    slots = parameter_slot_bindings(context.definition, named, (), process)
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
