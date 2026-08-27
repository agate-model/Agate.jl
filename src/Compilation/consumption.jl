function _consumption_resource_tracers(resources::Tuple, layout::ModelLayout)
    tracers = Symbol[]
    for resource in resources
        append!(tracers, getproperty(layout.component_tracers, resource))
    end
    return Tuple(tracers)
end

function _consumption_resources(
    ::PreferentialGrazing,
    resources::Tuple,
    layout::ModelLayout,
)
    return _realize_normalized_population_states(resources, layout)
end

function _consumption_resources(
    ::HeterotrophicConsumption,
    resources::Tuple,
    layout::ModelLayout,
)
    return _consumption_resource_tracers(resources, layout)
end

function _consumption_axis_positions(consumer_axis::Int, resource_axis::Int)
    return (
        consumer=_axis_position(consumer_axis),
        resource=_axis_position(resource_axis),
    )
end

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
        parameter_operand(slots.maximum_rate, context.plan, axis_positions),
        parameter_operand(slots.half_saturation, context.plan, axis_positions),
        parameter_operand(slots.palatability, context.plan, axis_positions),
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
        parameter_operand(slots.maximum_rate, context.plan, axis_positions),
        parameter_operand(slots.half_saturation, context.plan, axis_positions),
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
    consumer_tracers = _realize_normalized_population_states(
        named.facts.consumer_states, layout
    )
    resource_tracers = _consumption_resources(formulation, named.facts.resources, layout)
    slots = parameter_slot_bindings(context.definition, named, (), process)
    fluxes = Any[]

    for consumer_axis in eachindex(consumer_tracers)
        consumer = consumer_tracers[consumer_axis]
        for resource_axis in eachindex(resource_tracers)
            resource = resource_tracers[resource_axis]
            axis_positions = _consumption_axis_positions(consumer_axis, resource_axis)
            rate = _consumption_rate(
                formulation, slots, context, named, consumer, resource, axis_positions
            )
            assimilation = parameter_operand(slots.assimilation, context.plan, axis_positions)
            push!(
                fluxes,
                FluxSpec(resource, rate, Weight{-1}()),
                FluxSpec(consumer, rate, Weight{1}((assimilation,))),
            )
            if process.products isa Products
                append!(
                    fluxes,
                    _product_fluxes(
                        named, process.products, context, rate;
                        suffix=(RemainderOp((assimilation,)),),
                    ),
                )
            end
        end
    end
    return Tuple(fluxes)
end
