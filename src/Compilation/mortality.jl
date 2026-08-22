function _mortality_slots(
    definition::NormalizedModelDefinition, named::NamedProcess
)
    populations = named.process.populations
    context = length(populations) == 1 ? (population=only(populations),) : NamedTuple()
    return parameter_slot_bindings(
        definition, named, (), formulation(named.process); context
    )
end

function _mortality_rate(
    formulation,
    rate_binding::ParameterBinding,
    context::CommunityContext,
    population_axis::Int,
    population_index::Int,
    population_tracer::Symbol,
)
    axis_positions = (population=_axis_position(population_axis, population_index),)
    operands = (
        TracerOp{population_tracer}(),
        parameter_operand(rate_binding, context, axis_positions),
    )
    return RateElement(formulation, operands)
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Mortality}
    process = named.process
    population_tracers, population_indices = _realize_population_classes(
        named, process.populations, layout, context
    )
    slots = _mortality_slots(definition, named)
    fluxes = Any[]

    for population_axis in eachindex(population_tracers)
        population_index = population_indices[population_axis]
        rate = _mortality_rate(
            formulation(process), slots.rate, context, population_axis, population_index,
            population_tracers[population_axis],
        )
        push!(
            fluxes,
            FluxSpec(
                process_id(named), population_tracers[population_axis], rate, Weight{-1}()
            ),
        )
        if !isnothing(process.routing)
            append!(fluxes, _routing_fluxes(named, definition, process.routing, layout, rate))
        end
    end
    return Tuple(fluxes)
end
