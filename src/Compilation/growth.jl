function _growth_rate(
    named::CanonicalProcess,
    context::CompileContext,
    participant,
    scale_ref::Int,
)
    axis_positions = (plankton=participant.position,)
    rate_factors = _factor_ops(context, named, axis_positions)
    operands = (
        input_operand(context.layout, participant.tracer),
        parameter_operand(scale_ref, context, axis_positions),
    )
    return RateOp(formulation(named.process), operands; factors=rate_factors)
end

function _growth_resource_fluxes(
    named::CanonicalProcess,
    context::CompileContext,
    rate::RateOp,
    participant,
)
    process = named.process
    layout = context.layout
    fluxes = Any[
        FluxSpec(
            _scalar_component_target(layout, process.reference_resource),
            rate,
            Weight{-1}(),
        ),
    ]
    for (element, resource) in pairs(process.additional_resources)
        operands = _material_ratio_operand(context, participant.component, element, participant.position)
        isempty(operands) && continue
        push!(fluxes, FluxSpec(
            _scalar_component_target(layout, resource), rate, Weight{-1}(operands),
        ))
    end
    return Tuple(fluxes)
end

"""Derive biomass-gain and resource-loss fluxes for factorized growth."""
function process_fluxes(
    named::CanonicalProcess{Process}, context::CompileContext
) where {Process<:Growth}
    layout = context.layout
    participants = _realize_participants(named.semantic_facts.plankton_states, layout)
    scale_ref = named.binding_refs.process.maximum_rate
    fluxes = Any[]

    for participant in participants
        rate = _growth_rate(named, context, participant, scale_ref)
        product_targets = named.semantic_facts.product_targets
        if isnothing(product_targets)
            push!(fluxes, FluxSpec(participant.tracer, rate, Weight{1}()))
        else
            axis_positions = (plankton=participant.position,)
            product_fraction = parameter_operand(
                named.binding_refs.process.product_fraction, context, axis_positions
            )
            retained_fraction = ComplementOp((product_fraction,))
            push!(
                fluxes,
                FluxSpec(participant.tracer, rate, Weight{1}((retained_fraction,))),
            )
            for element in keys(first(values(product_targets)))
                ratio = _material_ratio_operand(
                    context, participant.component, element, participant.position
                )
                element !== named.semantic_facts.reference_element && isempty(ratio) && continue
                append!(fluxes, _product_fluxes_for_element(
                    named, product_targets, context, rate, element;
                    suffix=(product_fraction, ratio...),
                ))
            end
        end
        append!(fluxes, _growth_resource_fluxes(named, context, rate, participant))
    end
    return Tuple(fluxes)
end
