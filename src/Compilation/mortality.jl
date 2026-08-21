"""Realized mortality topology over concrete population classes and routing targets."""
struct MortalityTopology{TR,IX,R,E,D}
    population_tracers::TR
    population_indices::IX
    retained_target::R
    exported_target::E
    routed_targets::D
end

function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Mortality}
    process = named.process
    tracer_values, index_values = _realize_population_classes(
        named, process.populations, layout, context
    )

    routing = process.routing
    if isnothing(routing)
        return MortalityTopology(tracer_values, index_values, nothing, nothing, nothing)
    end
    if routing.formulation isa DirectRouting
        destination = _scalar_component_target(layout, routing.retained)
        return MortalityTopology(tracer_values, index_values, destination, nothing, nothing)
    elseif routing.formulation isa PartitionRouting
        retained = _scalar_component_target(layout, routing.retained)
        exported = _scalar_component_target(layout, routing.exported)
        return MortalityTopology(tracer_values, index_values, retained, exported, nothing)
    elseif routing.formulation isa DOMPOMRouting
        targets = _realize_dom_pom_routing(routing, layout)
        return MortalityTopology(tracer_values, index_values, nothing, nothing, targets)
    end
    throw(ArgumentError("unsupported mortality product routing $(typeof(routing.formulation))"))
end

function _mortality_slots(
    definition::NormalizedModelDefinition, named::NamedProcess
)
    populations = named.process.populations
    context = length(populations) == 1 ? (population=only(populations),) : NamedTuple()
    return parameter_slot_bindings(
        definition, named, (), formulation(named.process); context
    )
end

function _mortality_rate(formulation, rate_binding::ParameterBinding, population_index::Int)
    operands = (
        ClassOp{population_index}(),
        parameter_operand(rate_binding, population_index),
    )
    return RateElement(formulation, operands)
end

function _mortality_routed_fluxes(
    named::NamedProcess,
    definition::NormalizedModelDefinition,
    routed_targets::NamedTuple,
    rate::RateElement,
)
    routing = named.process.routing
    fraction = _routing_fraction_binding(definition, named, routing)
    fluxes = ()
    for route in (:DOM, :POM)
        targets = getproperty(routed_targets, route)
        for currency in keys(targets)
            target = getproperty(targets, currency)
            ratio = _routing_ratio_binding(definition, named, routing, currency)
            flux = FluxSpec(
                process_id(named),
                target,
                rate,
                _routing_weight(fraction, route; ratio),
            )
            fluxes = (fluxes..., flux)
        end
    end
    return fluxes
end

function process_fluxes(
    named::NamedProcess{P},
    topology::MortalityTopology,
    definition::NormalizedModelDefinition,
) where {P<:Mortality}
    slots = _mortality_slots(definition, named)
    form = formulation(named.process)
    routing = named.process.routing
    fraction = if isnothing(routing) || routing.formulation isa DirectRouting ||
                  !isnothing(topology.routed_targets)
        nothing
    else
        _routing_fraction_binding(definition, named, routing)
    end
    fluxes = ()

    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        rate = _mortality_rate(form, slots.rate, topology.population_indices[i])
        fluxes = (fluxes..., FluxSpec(process_id(named), tracer, rate, Weight{-1}()))

        if !isnothing(routing) && routing.formulation isa DirectRouting
            routed = FluxSpec(
                process_id(named), topology.retained_target, rate, Weight{1}()
            )
            fluxes = (fluxes..., routed)
        elseif !isnothing(topology.routed_targets)
            routed = _mortality_routed_fluxes(
                named, definition, topology.routed_targets, rate
            )
            fluxes = (fluxes..., routed...)
        elseif !isnothing(topology.retained_target)
            retained = FluxSpec(
                process_id(named),
                topology.retained_target,
                rate,
                _routing_weight(fraction, :retained),
            )
            exported = FluxSpec(
                process_id(named),
                topology.exported_target,
                rate,
                _routing_weight(fraction, :exported),
            )
            fluxes = (fluxes..., retained, exported)
        end
    end
    return fluxes
end

function process_fluxes(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Mortality}
    topology = realize_process_topology(named, layout, context)
    return process_fluxes(named, topology, definition)
end
