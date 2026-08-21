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

function _consumption_rate(
    formulation,
    slots,
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    layout::ComponentLayout,
    consumer_index::Int,
    consumer_axis::Int,
    resource::Symbol,
    resource_axis::Int,
    context::CommunityContext,
)
    axis_positions = (
        consumer=_axis_position(consumer_axis, consumer_index),
        resource=_axis_position(resource_axis),
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
    process.formulation isa HeterotrophicConsumption || throw(
        ArgumentError("unsupported consumption formulation $(typeof(process.formulation))"),
    )
    consumer_tracers, consumer_indices = _realize_population_classes(
        named, process.consumers, layout, context
    )
    resource_tracers = _consumption_resource_tracers(named, process.resources, layout)
    destination = isnothing(process.unassimilated_destination) ? nothing :
                  _scalar_component_target(layout, process.unassimilated_destination)
    slots = parameter_slot_bindings(definition, named, (), formulation(process))
    fluxes = ()

    for consumer_axis in eachindex(consumer_tracers)
        consumer = consumer_tracers[consumer_axis]
        consumer_index = consumer_indices[consumer_axis]
        for resource_axis in eachindex(resource_tracers)
            resource = resource_tracers[resource_axis]
            rate = _consumption_rate(
                formulation(process),
                slots,
                definition,
                named,
                layout,
                consumer_index,
                consumer_axis,
                resource,
                resource_axis,
                context,
            )
            axis_positions = (
                consumer=_axis_position(consumer_axis, consumer_index),
                resource=_axis_position(resource_axis),
            )
            assimilation = parameter_operand(slots.assimilation, context, axis_positions)
            fluxes = (
                fluxes...,
                FluxSpec(process_id(named), resource, rate, Weight{-1}()),
                FluxSpec(process_id(named), consumer, rate, Weight{1}((assimilation,))),
            )
            if !isnothing(destination)
                fluxes = (
                    fluxes...,
                    FluxSpec(
                        process_id(named),
                        destination,
                        rate,
                        Weight{1}((ComplementOp(assimilation),)),
                    ),
                )
            end
        end
    end
    return fluxes
end
