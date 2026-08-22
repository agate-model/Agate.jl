function _consumption_resource_tracers(
    named::NamedProcess, resources::Tuple, layout::ComponentLayout
)
    tracers = ()
    for resource in resources
        hasproperty(layout.component_tracers, resource) || throw(
            ArgumentError(
                "process :$(process_id(named)) references unrealized resource :$resource"
            ),
        )
        tracers = (tracers..., getproperty(layout.component_tracers, resource)...)
    end
    return tracers
end

function _consumption_resources(
    ::Union{IdealizedGrazing,PreferentialGrazing},
    named::NamedProcess,
    resources::Tuple,
    layout::ComponentLayout,
    context::CommunityContext,
)
    return _realize_population_classes(named, resources, layout, context)
end

function _consumption_resources(
    ::HeterotrophicConsumption,
    named::NamedProcess,
    resources::Tuple,
    layout::ComponentLayout,
    ::CommunityContext,
)
    return _consumption_resource_tracers(named, resources, layout), nothing
end

function _consumption_axis_positions(
    consumer_axis::Int,
    consumer_index::Int,
    resource_axis::Int,
    resource_index::Union{Nothing,Int},
)
    resource_position = isnothing(resource_index) ?
                        _axis_position(resource_axis) :
                        _axis_position(resource_axis, resource_index)
    return (
        consumer=_axis_position(consumer_axis, consumer_index),
        resource=resource_position,
    )
end

function _consumption_rate(
    formulation::Union{IdealizedGrazing,PreferentialGrazing},
    slots,
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    layout::ComponentLayout,
    consumer_index::Int,
    resource::Symbol,
    resource_index::Int,
    context::CommunityContext,
    axis_positions::NamedTuple,
)
    operands = (
        ClassOp{resource_index}(),
        ClassOp{consumer_index}(),
        parameter_operand(slots.maximum_rate, context, axis_positions),
        parameter_operand(slots.half_saturation, context, axis_positions),
        parameter_operand(slots.palatability, context, axis_positions),
    )
    rate_factors = _factor_elements(definition, named, layout, context, axis_positions)
    return RateElement(formulation, operands; factors=rate_factors)
end

function _consumption_rate(
    formulation::HeterotrophicConsumption,
    slots,
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    layout::ComponentLayout,
    consumer_index::Int,
    resource::Symbol,
    ::Nothing,
    context::CommunityContext,
    axis_positions::NamedTuple,
)
    operands = (
        TracerOp{resource}(),
        ClassOp{consumer_index}(),
        parameter_operand(slots.maximum_rate, context, axis_positions),
        parameter_operand(slots.half_saturation, context, axis_positions),
    )
    rate_factors = _factor_elements(definition, named, layout, context, axis_positions)
    return RateElement(formulation, operands; factors=rate_factors)
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Consumption}
    process = named.process
    formulation = process.formulation
    formulation isa Union{IdealizedGrazing,PreferentialGrazing,HeterotrophicConsumption} || throw(
        ArgumentError("unsupported consumption formulation $(typeof(formulation))"),
    )
    consumer_tracers, consumer_indices = _realize_population_classes(
        named, process.consumers, layout, context
    )
    resource_tracers, resource_indices = _consumption_resources(
        formulation, named, process.resources, layout, context
    )
    slots = parameter_slot_bindings(definition, named, (), formulation)
    fluxes = ()

    for consumer_axis in eachindex(consumer_tracers)
        consumer = consumer_tracers[consumer_axis]
        consumer_index = consumer_indices[consumer_axis]
        for resource_axis in eachindex(resource_tracers)
            resource = resource_tracers[resource_axis]
            resource_index = isnothing(resource_indices) ?
                             nothing : resource_indices[resource_axis]
            axis_positions = _consumption_axis_positions(
                consumer_axis, consumer_index, resource_axis, resource_index
            )
            rate = _consumption_rate(
                formulation,
                slots,
                definition,
                named,
                layout,
                consumer_index,
                resource,
                resource_index,
                context,
                axis_positions,
            )
            assimilation = parameter_operand(slots.assimilation, context, axis_positions)
            fluxes = (
                fluxes...,
                FluxSpec(process_id(named), resource, rate, Weight{-1}()),
                FluxSpec(process_id(named), consumer, rate, Weight{1}((assimilation,))),
            )
            isnothing(process.routing) || (
                fluxes = (
                    fluxes...,
                    _routing_fluxes(
                        named,
                        definition,
                        process.routing,
                        layout,
                        rate;
                        suffix=(ComplementOp(assimilation),),
                    )...,
                )
            )
        end
    end
    return fluxes
end
