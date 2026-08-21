function _grazing_rate(
    formulation,
    slots,
    consumer_index::Int,
    consumer_axis::Int,
    resource_index::Int,
    resource_axis::Int,
    context::CommunityContext,
)
    axis_positions = (
        consumer=_axis_position(consumer_axis, consumer_index),
        resource=_axis_position(resource_axis, resource_index),
    )
    operands = (
        ClassOp{resource_index}(),
        ClassOp{consumer_index}(),
        parameter_operand(slots.maximum_rate, context, axis_positions),
        parameter_operand(slots.half_saturation, context, axis_positions),
        parameter_operand(slots.palatability, context, axis_positions),
    )
    return RateElement(formulation, operands)
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Grazing}
    process = named.process
    process.formulation isa Union{IdealizedGrazing,PreferentialGrazing} || throw(
        ArgumentError("unsupported grazing formulation $(typeof(process.formulation))"),
    )
    isnothing(process.routing) || process.routing.formulation isa DOMPOMRouting || throw(
        ArgumentError(
            "unsupported grazing product routing $(typeof(process.routing.formulation))"
        ),
    )
    consumer_tracers, consumer_indices = _realize_population_classes(
        named, process.consumers, layout, context
    )
    resource_tracers, resource_indices = _realize_population_classes(
        named, process.resources, layout, context
    )
    destination = isnothing(process.unassimilated_destination) ? nothing :
                  _scalar_component_target(layout, process.unassimilated_destination)
    slots = parameter_slot_bindings(definition, named, (), formulation(process))
    fluxes = ()

    for consumer_axis in eachindex(consumer_tracers)
        consumer = consumer_tracers[consumer_axis]
        consumer_index = consumer_indices[consumer_axis]
        for resource_axis in eachindex(resource_tracers)
            resource = resource_tracers[resource_axis]
            rate = _grazing_rate(
                formulation(process),
                slots,
                consumer_index,
                consumer_axis,
                resource_indices[resource_axis],
                resource_axis,
                context,
            )
            axis_positions = (
                consumer=_axis_position(consumer_axis, consumer_index),
                resource=_axis_position(resource_axis, resource_indices[resource_axis]),
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
            elseif !isnothing(process.routing)
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
            end
        end
    end
    return fluxes
end
