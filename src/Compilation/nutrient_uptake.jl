function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:NutrientUptake}
    process = named.process
    facts = named.facts
    layout = context.layout
    target_tracers = state_tracers(layout, facts.target)
    reference_tracers = state_tracers(layout, facts.reference)
    resource = _scalar_component_target(layout, facts.resource)
    population_classes = component_classes(layout, process.population)
    slots = named.binding_refs.process
    fluxes = Any[]

    for population_axis in eachindex(target_tracers)
        axis_positions = (
            population=_axis_position(population_axis, population_classes[population_axis]),
        )
        operands = (
            input_operand(layout, resource),
            input_operand(layout, target_tracers[population_axis]),
            input_operand(layout, reference_tracers[population_axis]),
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
            FluxSpec(target_tracers[population_axis], rate, Weight{1}()),
        )
    end
    return Tuple(fluxes)
end
