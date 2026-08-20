"""Realized preferential-grazing topology over process-local interaction axes."""
struct GrazingTopology{CT,CI,RT,RI,D,R}
    consumer_tracers::CT
    consumer_indices::CI
    resource_tracers::RT
    resource_indices::RI
    unassimilated_target::D
    routing::R
end

function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Grazing}
    process = named.process
    process.formulation isa PreferentialGrazing || throw(
        ArgumentError("unsupported grazing formulation $(typeof(process.formulation))"),
    )

    consumer_tracers, consumer_indices = _realize_population_classes(
        named, process.consumers, layout, context
    )
    resource_tracers, resource_indices = _realize_population_classes(
        named, process.resources, layout, context
    )
    destination = if isnothing(process.unassimilated_destination)
        nothing
    else
        _scalar_component_target(layout, process.unassimilated_destination)
    end
    routing = if isnothing(process.routing)
        nothing
    elseif process.routing.formulation isa DOMPOMRouting
        _realize_dom_pom_routing(process.routing, layout)
    else
        throw(
            ArgumentError(
                "unsupported grazing product routing $(typeof(process.routing.formulation))"
            ),
        )
    end
    return GrazingTopology(
        consumer_tracers,
        consumer_indices,
        resource_tracers,
        resource_indices,
        destination,
        routing,
    )
end

function _grazing_rate(
    formulation,
    slots,
    consumer_index::Int,
    consumer_axis::Int,
    resource_index::Int,
    resource_axis::Int,
)
    operands = (
        ClassOp{resource_index}(),
        ClassOp{consumer_index}(),
        parameter_operand(slots.maximum_rate, consumer_index),
        parameter_operand(slots.half_saturation, consumer_index),
        parameter_operand(slots.palatability, consumer_axis, resource_axis),
    )
    return RateElement(formulation, operands)
end

function _grazing_routed_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    topology::GrazingTopology,
    rate::RateElement,
    assimilation,
)
    routing = named.process.routing
    isnothing(routing) && return ()
    fraction = _routing_fraction_binding(definition, named, routing)
    fluxes = ()
    for route in (:DOM, :POM)
        targets = getproperty(topology.routing, route)
        for currency in keys(targets)
            target = getproperty(targets, currency)
            ratio = _routing_ratio_binding(definition, named, routing, currency)
            weight = _routing_weight(
                fraction,
                route;
                ratio,
                suffix=(ComplementOp(assimilation),),
            )
            fluxes = (fluxes..., FluxSpec(process_id(named), target, rate, weight))
        end
    end
    return fluxes
end

function process_fluxes(
    named::NamedProcess{P},
    topology::GrazingTopology,
    definition::NormalizedModelDefinition,
) where {P<:Grazing}
    length(topology.consumer_tracers) == length(topology.consumer_indices) || throw(
        ArgumentError("grazing topology consumer tracer and index counts must match"),
    )
    length(topology.resource_tracers) == length(topology.resource_indices) || throw(
        ArgumentError("grazing topology resource tracer and index counts must match"),
    )

    slots = parameter_slot_bindings(definition, named, (), formulation(named.process))
    fluxes = ()
    for consumer_axis in eachindex(topology.consumer_tracers)
        consumer = topology.consumer_tracers[consumer_axis]
        consumer_index = topology.consumer_indices[consumer_axis]
        for resource_axis in eachindex(topology.resource_tracers)
            resource = topology.resource_tracers[resource_axis]
            resource_index = topology.resource_indices[resource_axis]
            rate = _grazing_rate(
                formulation(named.process),
                slots,
                consumer_index,
                consumer_axis,
                resource_index,
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
            routed = _grazing_routed_fluxes(
                named, definition, topology, rate, assimilation
            )
            fluxes = (fluxes..., routed...)
        end
    end
    return fluxes
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Grazing}
    topology = realize_process_topology(named, layout, context)
    return process_fluxes(named, topology, definition)
end
