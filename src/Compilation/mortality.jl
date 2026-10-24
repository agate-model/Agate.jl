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
    named::CanonicalProcess{Process}, context::CompileContext
) where {Process<:Mortality}
    process = named.process
    layout = context.layout
    fluxes = Any[]

    for reference in named.semantic_facts.plankton_states
        slots = getproperty(named.binding_refs.process, reference.plankton)
        state_refs = getproperty(named.semantic_facts.state_sets, reference.plankton)
        materials = getproperty(named.semantic_facts.material_sources, reference.plankton)
        for participant in _realize_participants((reference,), layout)
            for state_ref in state_refs
                tracer = state_tracer(layout, state_ref, participant.component_index)
                rate = _mortality_rate(
                    formulation(process), slots.rate, context, participant, tracer
                )
                push!(fluxes, FluxSpec(tracer, rate, Weight{-1}()))
                isnothing(named.semantic_facts.product_targets) && continue
                for (element, source) in pairs(materials)
                    source.state === state_ref.state || continue
                    ratio = _material_ratio_operand(
                        context, participant.component, element, participant.position
                    )
                    append!(fluxes, _product_fluxes_for_element(
                        named, named.semantic_facts.product_targets, context, rate, element;
                        suffix=ratio,
                    ))
                end
            end
        end
    end
    return Tuple(fluxes)
end
