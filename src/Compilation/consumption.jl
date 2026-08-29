function _consumption_rate(
    formulation::PreferentialGrazing,
    slots,
    context::CompileContext,
    named::NamedProcess,
    inventory::Symbol,
    reference_resource::Symbol,
    consumer::Symbol,
    axis_positions::NamedTuple,
)
    operands = (
        input_operand(context.layout, inventory),
        input_operand(context.layout, reference_resource),
        input_operand(context.layout, consumer),
        parameter_operand(slots.maximum_rate, context, axis_positions),
        parameter_operand(slots.half_saturation, context, axis_positions),
        parameter_operand(slots.palatability, context, axis_positions),
    )
    rate_factors = _factor_ops(context, named, axis_positions)
    return RateOp(formulation, operands; factors=rate_factors)
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
    rate_factors = _factor_ops(context, named, axis_positions)
    return RateOp(formulation, operands; factors=rate_factors)
end

function _append_unassimilated_products!(
    fluxes,
    named,
    context,
    rate,
    assimilation,
    element,
)
    isnothing(named.facts.product_targets) && return nothing
    suffix = (RemainderOp((assimilation,)),)
    if named.facts.product_mode === :state
        append!(
            fluxes,
            _product_fluxes_for_element(
                named, named.facts.product_targets, context, rate, element; suffix
            ),
        )
    else
        append!(
            fluxes,
            _product_fluxes(named, named.facts.product_targets, context, rate; suffix),
        )
    end
    return nothing
end

function _living_consumption_fluxes!(
    fluxes,
    named::NamedProcess,
    context::CompileContext,
    consumer,
    resource,
    slots,
    axis_positions,
)
    layout = context.layout
    state_refs = getproperty(named.facts.resource_state_sets, resource.component)
    state_elements = getproperty(named.facts.resource_state_elements, resource.component)
    consumer_element_states = getproperty(
        named.facts.consumer_element_states, consumer.component
    )
    assimilation = parameter_operand(slots.assimilation, context, axis_positions)

    for state_ref in state_refs
        resource_tracer = state_tracer(layout, state_ref, resource.component_index)
        rate = _consumption_rate(
            named.process.formulation,
            slots,
            context,
            named,
            resource_tracer,
            resource.tracer,
            consumer.tracer,
            axis_positions,
        )
        push!(fluxes, FluxSpec(resource_tracer, rate, Weight{-1}()))

        state_element_value = getproperty(state_elements, state_ref.state)
        isnothing(state_element_value) && continue
        consumer_ref = getproperty(consumer_element_states, state_element_value)
        consumer_tracer = state_tracer(layout, consumer_ref, consumer.component_index)
        push!(fluxes, FluxSpec(consumer_tracer, rate, Weight{1}((assimilation,))))

        if named.facts.product_mode === :stoichiometric &&
           state_element_value !== named.facts.reference_element
            continue
        end
        _append_unassimilated_products!(
            fluxes, named, context, rate, assimilation, state_element_value
        )
    end
    return nothing
end

function _heterotrophic_consumption_fluxes!(
    fluxes,
    named::NamedProcess,
    context::CompileContext,
    consumer,
    resource,
    slots,
    axis_positions,
)
    layout = context.layout
    rate = _consumption_rate(
        named.process.formulation,
        slots,
        context,
        named,
        consumer.tracer,
        resource.tracer,
        axis_positions,
    )
    assimilation = parameter_operand(slots.assimilation, context, axis_positions)
    consumer_element_states = getproperty(
        named.facts.consumer_element_states, consumer.component
    )
    consumer_ref = getproperty(consumer_element_states, named.facts.reference_element)
    consumer_tracer = state_tracer(layout, consumer_ref, consumer.component_index)
    push!(
        fluxes,
        FluxSpec(resource.tracer, rate, Weight{-1}()),
        FluxSpec(consumer_tracer, rate, Weight{1}((assimilation,))),
    )
    _append_unassimilated_products!(
        fluxes,
        named,
        context,
        rate,
        assimilation,
        named.facts.reference_element,
    )
    return nothing
end

function process_fluxes(
    named::NamedProcess{P}, context::CompileContext
) where {P<:Consumption}
    form = named.process.formulation
    layout = context.layout
    consumers = _realize_participants(named.facts.consumer_states, layout)
    resources = _realize_participants(named.facts.resources, layout)
    slots = named.binding_refs.process
    fluxes = Any[]

    for consumer in consumers, resource in resources
        axis_positions = (consumer=consumer.position, resource=resource.position)
        if form isa PreferentialGrazing
            _living_consumption_fluxes!(
                fluxes, named, context, consumer, resource, slots, axis_positions
            )
        else
            _heterotrophic_consumption_fluxes!(
                fluxes, named, context, consumer, resource, slots, axis_positions
            )
        end
    end
    return Tuple(fluxes)
end
