function _mortality_slots(
    context::CompileContext, named::NamedProcess, population::Symbol
)
    return parameter_slot_bindings(
        context.definition, named, (), named.process; context=(population=population,)
    )
end

function _mortality_rate(
    formulation,
    rate_binding::ParameterBinding,
    context::CompileContext,
    population_axis::Int,
    population_class::Symbol,
    population_tracer::Symbol,
)
    axis_positions = (population=_axis_position(population_axis, population_class),)
    operands = (
        input_operand(context.layout, population_tracer),
        parameter_operand(rate_binding, context.plan, axis_positions),
    )
    return RateElement(formulation, operands)
end

function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:Mortality}
    process = named.process
    fluxes = Any[]

    for reference in named.facts.population_states
        population = reference.population
        population_tracers = _realize_population_states((reference,), context.layout)
        population_classes = component_classes(context.layout, population)
        slots = _mortality_slots(context, named, population)
        for population_axis in eachindex(population_tracers)
            rate = _mortality_rate(
                formulation(process), slots.rate, context, population_axis,
                population_classes[population_axis], population_tracers[population_axis],
            )
            push!(
                fluxes,
                FluxSpec(population_tracers[population_axis], rate, Weight{-1}()),
            )
            if !isnothing(named.facts.product_targets)
                append!(
                    fluxes, _product_fluxes(named, named.facts.product_targets, context, rate)
                )
            end
        end
    end
    return Tuple(fluxes)
end
