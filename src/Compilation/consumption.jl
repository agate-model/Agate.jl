function _consumption_rate(
    formulation::PreferentialGrazing,
    slots,
    context::CompileContext,
    named::CanonicalProcess,
    inventory::Symbol,
    reference_resource::Symbol,
    consumer::Symbol,
    axis_positions::NamedTuple,
    shared_operands::Tuple,
)
    operands = (
        input_operand(context.layout, inventory),
        input_operand(context.layout, reference_resource),
        input_operand(context.layout, consumer),
        parameter_operand(slots.maximum_rate, context, axis_positions),
        parameter_operand(slots.half_saturation, context, axis_positions),
        parameter_operand(slots.palatability, context, axis_positions),
        shared_operands...,
    )
    rate_factors = _factor_ops(context, named, axis_positions)
    return RateOp(formulation, operands; factors=rate_factors)
end

function _consumption_rate(
    formulation::HeterotrophicConsumption,
    slots,
    context::CompileContext,
    named::CanonicalProcess,
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
    isnothing(named.semantic_facts.product_targets) && return nothing
    suffix = (ComplementOp((assimilation,)),)
    if named.semantic_facts.product_mode === :state
        append!(
            fluxes,
            _product_fluxes_for_element(
                named, named.semantic_facts.product_targets, context, rate, element; suffix
            ),
        )
    else
        append!(
            fluxes,
            _product_fluxes(named, named.semantic_facts.product_targets, context, rate; suffix),
        )
    end
    return nothing
end

function _living_consumption_fluxes!(
    fluxes,
    named::CanonicalProcess,
    context::CompileContext,
    consumer,
    resource,
    slots,
    axis_positions,
    shared_operands::Tuple,
)
    layout = context.layout
    state_refs = getproperty(named.semantic_facts.resource_state_sets, resource.component)
    state_elements = getproperty(named.semantic_facts.resource_state_elements, resource.component)
    consumer_element_states = getproperty(
        named.semantic_facts.consumer_element_states, consumer.component
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
            shared_operands,
        )
        push!(fluxes, FluxSpec(resource_tracer, rate, Weight{-1}()))

        state_element_value = getproperty(state_elements, state_ref.state)
        isnothing(state_element_value) && continue
        consumer_ref = getproperty(consumer_element_states, state_element_value)
        consumer_tracer = state_tracer(layout, consumer_ref, consumer.component_index)
        push!(fluxes, FluxSpec(consumer_tracer, rate, Weight{1}((assimilation,))))

        if named.semantic_facts.product_mode === :stoichiometric &&
           state_element_value !== named.semantic_facts.reference_element
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
    named::CanonicalProcess,
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
        named.semantic_facts.consumer_element_states, consumer.component
    )
    consumer_ref = getproperty(consumer_element_states, named.semantic_facts.reference_element)
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
        named.semantic_facts.reference_element,
    )
    return nothing
end

function process_fluxes(
    named::CanonicalProcess{Process}, context::CompileContext
) where {Process<:Consumption}
    form = named.process.formulation
    layout = context.layout
    consumers = _realize_participants(named.semantic_facts.consumer_states, layout)
    resources = _realize_participants(named.semantic_facts.resources, layout)
    slots = named.binding_refs.process
    fluxes = Any[]

    if form isa PreferentialGrazing
        for consumer in consumers
            reference_resources = Tuple(
                input_operand(layout, resource.tracer) for resource in resources
            )
            palatabilities = Tuple(
                parameter_operand(
                    slots.palatability,
                    context,
                    (consumer=consumer.position, resource=resource.position),
                ) for resource in resources
            )
            # Keep consumer-level prey reductions as scalar IR nodes. Materializing the full
            # evaluated prey/palatability tuples in every edge rate causes generated-code
            # growth to become pathological for richer food webs.
            palatable_biomass = WeightedPowerSumOp{1}(reference_resources, palatabilities)
            switching_exponent = form.switching_exponent
            shared_operands = if switching_exponent == 1
                (palatable_biomass,)
            else
                switching_weights = WeightedPowerSumOp{switching_exponent}(
                    reference_resources, palatabilities
                )
                (palatable_biomass, switching_weights)
            end

            for resource in resources
                axis_positions = (consumer=consumer.position, resource=resource.position)
                _living_consumption_fluxes!(
                    fluxes,
                    named,
                    context,
                    consumer,
                    resource,
                    slots,
                    axis_positions,
                    shared_operands,
                )
            end
        end
    else
        for consumer in consumers, resource in resources
            axis_positions = (consumer=consumer.position, resource=resource.position)
            _heterotrophic_consumption_fluxes!(
                fluxes, named, context, consumer, resource, slots, axis_positions
            )
        end
    end
    return Tuple(fluxes)
end
