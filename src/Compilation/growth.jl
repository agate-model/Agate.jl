function _growth_rate(
    named::NamedProcess,
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
    named::NamedProcess,
    context::CompileContext,
    rate::RateOp,
)
    facts = named.facts
    layout = context.layout
    fluxes = Any[
        FluxSpec(
            _scalar_component_target(layout, facts.reference_source),
            rate,
            Weight{-1}(),
        ),
    ]
    for (element, resource) in pairs(facts.additional_resources)
        ratio_ref = getproperty(named.binding_refs.stoichiometry, element).ratio
        push!(
            fluxes,
            FluxSpec(
                _scalar_component_target(layout, resource),
                rate,
                Weight{-1}((parameter_operand(ratio_ref, context),)),
            ),
        )
    end
    return Tuple(fluxes)
end

"""Derive biomass-gain and resource-loss fluxes for factorized growth."""
function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:Growth}
    layout = context.layout
    participants = _realize_participants(named.facts.plankton_states, layout)
    scale_ref = named.binding_refs.process.maximum_rate
    fluxes = Any[]

    for participant in participants
        rate = _growth_rate(named, context, participant, scale_ref)
        push!(fluxes, FluxSpec(participant.tracer, rate, Weight{1}()))
        append!(fluxes, _growth_resource_fluxes(named, context, rate))
    end
    return Tuple(fluxes)
end
