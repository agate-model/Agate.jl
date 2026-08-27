function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:NutrientUptake}
    process = named.process
    facts = named.facts
    layout = context.layout
    participants = _realize_participants((facts.target,), layout)
    reference_tracers = state_tracers(layout, facts.reference)
    resource = _scalar_component_target(layout, facts.resource)
    slots = named.binding_refs.process
    fluxes = Any[]

    for participant in participants
        axis_positions = (population=participant.position,)
        operands = (
            input_operand(layout, resource),
            input_operand(layout, participant.tracer),
            input_operand(layout, reference_tracers[participant.position.local_index]),
            parameter_operand(slots.maximum_rate, context, axis_positions),
            parameter_operand(slots.K, context, axis_positions),
            parameter_operand(slots.minimum_quota, context, axis_positions),
            parameter_operand(slots.maximum_quota, context, axis_positions),
            parameter_operand(slots.hill, context, axis_positions),
        )
        rate = RateElement(formulation(process), operands)
        push!(
            fluxes,
            FluxSpec(resource, rate, Weight{-1}()),
            FluxSpec(participant.tracer, rate, Weight{1}()),
        )
    end
    return Tuple(fluxes)
end
