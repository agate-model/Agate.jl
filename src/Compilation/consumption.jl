function _consumption_rate(
    formulation::PreferentialGrazing,
    slots,
    context::CompileContext,
    named::NamedProcess,
    consumer::Symbol,
    resource::Symbol,
    axis_positions::NamedTuple,
)
    operands = (
        input_operand(context.layout, resource),
        input_operand(context.layout, consumer),
        parameter_operand(slots.maximum_rate, context, axis_positions),
        parameter_operand(slots.half_saturation, context, axis_positions),
        parameter_operand(slots.palatability, context, axis_positions),
    )
    rate_factors = _factor_elements(context, named, axis_positions)
    return RateElement(formulation, operands; factors=rate_factors)
end

function _consumption_rate(
    formulation::HeterotrophicConsumption,
    slots,
    context::CompileContext,
    named::NamedProcess,
    consumer::Symbol,
    resource::Symbol,
    axis_positions::NamedTuple,
)
    operands = (
        input_operand(context.layout, resource),
        input_operand(context.layout, consumer),
        parameter_operand(slots.maximum_rate, context, axis_positions),
        parameter_operand(slots.half_saturation, context, axis_positions),
    )
    rate_factors = _factor_elements(context, named, axis_positions)
    return RateElement(formulation, operands; factors=rate_factors)
end

function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:Consumption}
    process = named.process
    formulation = process.formulation
    layout = context.layout
    consumers = _realize_participants(named.facts.consumer_states, layout)
    resources = _realize_participants(named.facts.resources, layout)
    slots = named.binding_refs.process
    fluxes = Any[]

    for consumer in consumers
        for resource in resources
            axis_positions = (consumer=consumer.position, resource=resource.position)
            rate = _consumption_rate(
                formulation, slots, context, named, consumer.tracer, resource.tracer,
                axis_positions,
            )
            assimilation = parameter_operand(slots.assimilation, context, axis_positions)
            push!(
                fluxes,
                FluxSpec(resource.tracer, rate, Weight{-1}()),
                FluxSpec(consumer.tracer, rate, Weight{1}((assimilation,))),
            )
            if !isnothing(named.facts.product_targets)
                append!(
                    fluxes,
                    _product_fluxes(
                        named, named.facts.product_targets, context, rate;
                        suffix=(RemainderOp((assimilation,)),),
                    ),
                )
            end
        end
    end
    return Tuple(fluxes)
end
