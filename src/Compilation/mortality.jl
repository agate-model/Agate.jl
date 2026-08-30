function _mortality_rate(
    formulation,
    rate_ref::Int,
    context::CompileContext,
    participant,
    tracer::Symbol,
)
    axis_positions = (plankton=participant.position,)
    operands = (
        input_operand(context.layout, tracer),
        input_operand(context.layout, participant.tracer),
        parameter_operand(rate_ref, context, axis_positions),
    )
    return RateOp(formulation, operands)
end

function process_fluxes(
    named::CanonicalProcess{P}, context::CompileContext
) where {P<:Mortality}
    process = named.process
    layout = context.layout
    fluxes = Any[]

    for reference in named.semantic_facts.plankton_states
        slots = getproperty(named.binding_refs.process, reference.plankton)
        state_refs = getproperty(named.semantic_facts.state_sets, reference.plankton)
        state_elements = getproperty(named.semantic_facts.state_elements, reference.plankton)
        for participant in _realize_participants((reference,), layout)
            for state_ref in state_refs
                tracer = state_tracer(layout, state_ref, participant.component_index)
                rate = _mortality_rate(
                    formulation(process), slots.rate, context, participant, tracer
                )
                push!(fluxes, FluxSpec(tracer, rate, Weight{-1}()))

                isnothing(named.semantic_facts.product_targets) && continue
                state_element_value = getproperty(state_elements, state_ref.state)
                isnothing(state_element_value) && continue
                if named.semantic_facts.product_mode === :state
                    append!(
                        fluxes,
                        _product_fluxes_for_element(
                            named,
                            named.semantic_facts.product_targets,
                            context,
                            rate,
                            state_element_value,
                        ),
                    )
                elseif state_element_value === _state_element_for_reference(named, reference)
                    append!(
                        fluxes,
                        _product_fluxes(named, named.semantic_facts.product_targets, context, rate),
                    )
                end
            end
        end
    end
    return Tuple(fluxes)
end

function _state_element_for_reference(named::CanonicalProcess, reference::PlanktonStateRef)
    state_elements = getproperty(named.semantic_facts.state_elements, reference.plankton)
    return getproperty(state_elements, reference.state)
end
