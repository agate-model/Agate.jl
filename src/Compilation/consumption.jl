"""Realized consumer-by-resource topology for heterotrophic consumption."""
struct ConsumptionTopology{CT,CI,RT,D,L}
    consumer_tracers::CT
    consumer_indices::CI
    resource_tracers::RT
    unassimilated_target::D
    layout::L
end

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

function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Consumption}
    process = named.process
    process.formulation isa HeterotrophicConsumption || throw(
        ArgumentError("unsupported consumption formulation $(typeof(process.formulation))"),
    )
    consumer_tracers, consumer_indices = _realize_population_classes(
        named, process.consumers, layout, context
    )
    resources = _consumption_resource_tracers(named, process.resources, layout)
    destination = isnothing(process.unassimilated_destination) ? nothing :
                  _scalar_component_target(layout, process.unassimilated_destination)
    return ConsumptionTopology(consumer_tracers, consumer_indices, resources, destination, layout)
end

function _consumption_rate(
    formulation,
    slots,
    definition::NormalizedModelDefinition,
    named::NamedProcess,
    topology::ConsumptionTopology,
    consumer_index::Int,
    consumer_axis::Int,
    resource::Symbol,
    resource_axis::Int,
)
    operands = (
        TracerOp{resource}(),
        ClassOp{consumer_index}(),
        parameter_operand(slots.maximum_rate, consumer_axis),
        parameter_operand(slots.half_saturation, resource_axis),
    )
    axis_indices = (consumer=consumer_axis, resource=resource_axis)
    rate_factors = _factor_elements(definition, named, topology.layout, axis_indices)
    return RateElement(formulation, operands; factors=rate_factors)
end

function process_fluxes(
    named::NamedProcess{P},
    topology::ConsumptionTopology,
    definition::NormalizedModelDefinition,
) where {P<:Consumption}
    slots = parameter_slot_bindings(definition, named, (), formulation(named.process))
    fluxes = ()
    for consumer_axis in eachindex(topology.consumer_tracers)
        consumer = topology.consumer_tracers[consumer_axis]
        consumer_index = topology.consumer_indices[consumer_axis]
        for resource_axis in eachindex(topology.resource_tracers)
            resource = topology.resource_tracers[resource_axis]
            rate = _consumption_rate(
                formulation(named.process),
                slots,
                definition,
                named,
                topology,
                consumer_index,
                consumer_axis,
                resource,
                resource_axis,
            )
            assimilation = parameter_operand(
                slots.assimilation, consumer_axis, resource_axis
            )
            loss = FluxSpec(process_id(named), resource, rate, Weight{-1}())
            gain = FluxSpec(process_id(named), consumer, rate, Weight{1}((assimilation,)))
            fluxes = (fluxes..., loss, gain)
            if !isnothing(topology.unassimilated_target)
                unassimilated = FluxSpec(
                    process_id(named),
                    topology.unassimilated_target,
                    rate,
                    Weight{1}((ComplementOp(assimilation),)),
                )
                fluxes = (fluxes..., unassimilated)
            end
        end
    end
    return fluxes
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Consumption}
    topology = realize_process_topology(named, layout, context)
    return process_fluxes(named, topology, definition)
end
