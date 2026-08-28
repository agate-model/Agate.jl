function _mortality_rate(
    formulation,
    rate_ref::Int,
    context::CompileContext,
    participant,
)
    axis_positions = (population=participant.position,)
    operands = (
        input_operand(context.layout, participant.tracer),
        parameter_operand(rate_ref, context, axis_positions),
    )
    return RateElement(formulation, operands)
end

function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:Mortality}
    process = named.process
    fluxes = Any[]

    for reference in named.facts.population_states
        slots = getproperty(named.binding_refs.process, reference.population)
        for participant in _realize_participants((reference,), context.layout)
            rate = _mortality_rate(formulation(process), slots.rate, context, participant)
            push!(fluxes, FluxSpec(participant.tracer, rate, Weight{-1}()))
            if !isnothing(named.facts.product_targets)
                append!(
                    fluxes, _product_fluxes(named, named.facts.product_targets, context, rate)
                )
            end
        end
    end
    return Tuple(fluxes)
end
