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
    inventory_operand = input_operand(context.layout, inventory)
    common_operands = (
        input_operand(context.layout, consumer),
        parameter_operand(slots.maximum_rate, context, axis_positions),
        parameter_operand(slots.half_saturation, context, axis_positions),
        parameter_operand(slots.palatability, context, axis_positions),
        shared_operands...,
    )
    reference_operands = if formulation.switching_exponent == 1
        ()
    else
        (input_operand(context.layout, reference_resource),)
    end
    operands = (inventory_operand, reference_operands..., common_operands...)
    rate_factors = _factor_ops(context, named, axis_positions)
    return RateOp(formulation, operands; factors=rate_factors)
end

function _consumption_rate(
    formulation::LinearGrazing,
    slots,
    context::CompileContext,
    named::CanonicalProcess,
    inventory::Symbol,
    _reference_resource::Symbol,
    consumer::Symbol,
    axis_positions::NamedTuple,
    _shared_operands::Tuple,
)
    operands = (
        input_operand(context.layout, inventory),
        input_operand(context.layout, consumer),
        parameter_operand(slots.rate, context, axis_positions),
        parameter_operand(slots.palatability, context, axis_positions),
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
    shared_operands::Tuple,
)
    operands = (
        input_operand(context.layout, resource),
        input_operand(context.layout, consumer),
        parameter_operand(slots.maximum_rate, context, axis_positions),
        parameter_operand(slots.half_saturation, context, axis_positions),
        parameter_operand(slots.substrate_preference, context, axis_positions),
        shared_operands...,
    )
    rate_factors = _factor_ops(context, named, axis_positions)
    return RateOp(formulation, operands; factors=rate_factors)
end

function _append_unassimilated_products!(
    fluxes, named, context, rate, assimilation, element; ratio=(),
)
    isnothing(named.semantic_facts.product_targets) && return nothing
    append!(fluxes, _product_fluxes_for_element(
        named, named.semantic_facts.product_targets, context, rate, element;
        suffix=(ComplementOp((assimilation,)), ratio...),
    ))
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
    source_materials = getproperty(named.semantic_facts.resource_element_states, resource.component)
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

        for (element, source) in pairs(source_materials)
            source.state === state_ref.state || continue
            ratio = _material_ratio_operand(
                context, resource.component, element, resource.position
            )
            consumer_inventory = getproperty(consumer_element_states, element)
            consumer_tracer = state_tracer(
                layout, consumer_inventory, consumer.component_index
            )
            push!(fluxes, FluxSpec(
                consumer_tracer, rate, Weight{1}((assimilation, ratio...)),
            ))
            _append_unassimilated_products!(
                fluxes, named, context, rate, assimilation, element; ratio,
            )
        end
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
    shared_operands::Tuple,
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
        shared_operands,
    )
    assimilation = parameter_operand(slots.assimilation, context, axis_positions)
    consumer_element_states = getproperty(
        named.semantic_facts.consumer_element_states, consumer.component
    )
    consumer_inventory = getproperty(consumer_element_states, named.semantic_facts.reference_element)
    consumer_tracer = state_tracer(layout, consumer_inventory, consumer.component_index)
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

function _shared_consumption_operands(
    formulation::PreferentialGrazing, consumer, resources, slots, context::CompileContext
)
    layout = context.layout
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
    total_palatable_biomass = ProductPowerSumOp{1}(reference_resources, palatabilities)
    formulation.switching_exponent == 1 && return (total_palatable_biomass,)
    switching_weight_sum = ProductPowerSumOp{formulation.switching_exponent}(
        reference_resources, palatabilities
    )
    return (total_palatable_biomass, switching_weight_sum)
end

_shared_consumption_operands(
    ::LinearGrazing, consumer, resources, slots, context::CompileContext
) = ()

function _shared_consumption_operands(
    ::HeterotrophicConsumption, consumer, resources, slots, context::CompileContext
)
    layout = context.layout
    resource_operands = Tuple(
        input_operand(layout, resource.tracer) for resource in resources
    )
    substrate_affinities = Tuple(begin
        axis_positions = (consumer=consumer.position, resource=resource.position)
        QuotientOp(
            parameter_operand(slots.substrate_preference, context, axis_positions),
            parameter_operand(slots.half_saturation, context, axis_positions),
        )
    end for resource in resources)
    total_substrate_availability = ProductPowerSumOp{1}(
        resource_operands, substrate_affinities
    )
    return (total_substrate_availability,)
end

_consumption_flux_appender(::PreferentialGrazing) = _living_consumption_fluxes!
_consumption_flux_appender(::LinearGrazing) = _living_consumption_fluxes!
_consumption_flux_appender(::HeterotrophicConsumption) = _heterotrophic_consumption_fluxes!

function process_fluxes(
    named::CanonicalProcess{Process}, context::CompileContext
) where {Process<:Consumption}
    layout = context.layout
    formulation = named.process.formulation
    # A1.1 validates realized implicit ratios; A1.2 will lower their conservative transfers.
    if uses_living_interactions(formulation) && any(
        consumer -> !isempty(getproperty(context.definition.composition_refs, consumer)),
        named.process.consumers,
    )
        throw(ArgumentError(
            "process :$(process_id(named)) requires balancing transfers for " *
            "implicit consumer composition during grazing",
        ))
    end
    consumers = _realize_participants(named.semantic_facts.consumer_states, layout)
    resources = _realize_participants(named.semantic_facts.resources, layout)
    slots = named.binding_refs.process
    append_fluxes! = _consumption_flux_appender(formulation)
    fluxes = Any[]

    for consumer in consumers
        shared_operands = _shared_consumption_operands(
            formulation, consumer, resources, slots, context
        )
        for resource in resources
            axis_positions = (consumer=consumer.position, resource=resource.position)
            append_fluxes!(
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

    return Tuple(fluxes)
end
