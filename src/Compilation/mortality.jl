"""Resolved parameter names needed to compile one mortality process instance."""
struct MortalityParameterBinding{R}
    rate::Symbol
    routing_fraction::Union{Nothing,Symbol}
    routing::R
end

MortalityParameterBinding(
    rate::Symbol; routing_fraction::Union{Nothing,Symbol}=nothing
) = MortalityParameterBinding(rate, routing_fraction, nothing)

function MortalityParameterBinding(
    definition::NormalizedModelDefinition, id::Symbol
)
    hasproperty(definition.processes, id) || throw(
        ArgumentError("normalized model has no process :$id"),
    )
    named = getproperty(definition.processes, id)
    named.process isa Mortality || throw(
        ArgumentError("process :$id is not a Mortality process"),
    )

    rate = parameter_name(definition, _parameter_requirement(named, (), :rate))
    routing = named.process.routing
    isnothing(routing) && return MortalityParameterBinding(rate)
    if routing.formulation isa PartitionRouting
        routing_fraction = parameter_name(
            definition, _parameter_requirement(named, (:routing,), :export_fraction)
        )
        return MortalityParameterBinding(rate; routing_fraction)
    elseif routing.formulation isa DOMPOMRouting
        organic = _dom_pom_routing_binding(definition, named, routing)
        return MortalityParameterBinding(rate, organic.fraction, organic)
    end
    throw(ArgumentError("unsupported mortality product routing $(typeof(routing.formulation))"))
end

"""Realized mortality topology over concrete population classes and routing targets."""
struct MortalityTopology{TR,IX,R,E}
    population_tracers::TR
    population_indices::IX
    retained_target::R
    exported_target::E
end

struct DOMPOMMortalityTopology{TR,IX,R}
    population_tracers::TR
    population_indices::IX
    routing::R
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
        return MortalityTopology(tracer_values, index_values, nothing, nothing)
    end
    if routing.formulation isa PartitionRouting
        retained = _scalar_component_target(layout, routing.retained)
        exported = _scalar_component_target(layout, routing.exported)
        return MortalityTopology(tracer_values, index_values, retained, exported)
    elseif routing.formulation isa DOMPOMRouting
        return DOMPOMMortalityTopology(
            tracer_values, index_values, _realize_dom_pom_routing(routing, layout)
        )
    end
    throw(ArgumentError("unsupported mortality product routing $(typeof(routing.formulation))"))
end

function _validate_binding(process::Mortality, binding::MortalityParameterBinding)
    isnothing(process.routing) && return nothing
    isnothing(binding.routing_fraction) && throw(
        ArgumentError("routed mortality requires a routing-fraction parameter binding"),
    )
    return nothing
end

function _mortality_rate(formulation, population_index::Int, parameter::Symbol)
    operands = (ClassOp{population_index}(), VecParamOp{parameter,population_index}())
    return RateElement(formulation, operands)
end

function _mortality_routed_fluxes(
    named::NamedProcess,
    topology::DOMPOMMortalityTopology,
    binding::MortalityParameterBinding,
    rate::RateElement,
)
    routing_binding = binding.routing
    isnothing(routing_binding) && throw(
        ArgumentError("DOM/POM-routed mortality requires a routing parameter binding"),
    )
    fluxes = ()
    for route in (:DOM, :POM)
        targets = getproperty(topology.routing, route)
        for currency in keys(targets)
            target = getproperty(targets, currency)
            ratio = _routing_ratio_parameter(topology.routing, routing_binding, currency)
            flux = FluxSpec(
                process_id(named),
                target,
                rate,
                _routing_weight(routing_binding.fraction, route; ratio),
            )
            fluxes = (fluxes..., flux)
        end
    end
    return fluxes
end

function process_fluxes(
    named::NamedProcess{P}, topology::MortalityTopology, binding::MortalityParameterBinding
) where {P<:Mortality}
    length(topology.population_tracers) == length(topology.population_indices) || throw(
        ArgumentError("mortality topology tracer and index counts must match"),
    )
    _validate_binding(named.process, binding)
    form = formulation(named.process)
    fluxes = ()

    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        rate = _mortality_rate(form, topology.population_indices[i], binding.rate)
        fluxes = (fluxes..., FluxSpec(process_id(named), tracer, rate, Weight{-1}()))

        if !isnothing(topology.retained_target)
            fraction = something(binding.routing_fraction)
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
    topology::DOMPOMMortalityTopology,
    binding::MortalityParameterBinding,
) where {P<:Mortality}
    length(topology.population_tracers) == length(topology.population_indices) || throw(
        ArgumentError("mortality topology tracer and index counts must match"),
    )
    fluxes = ()
    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        rate = _mortality_rate(
            formulation(named.process), topology.population_indices[i], binding.rate
        )
        routed = _mortality_routed_fluxes(named, topology, binding, rate)
        fluxes = (
            fluxes...,
            FluxSpec(process_id(named), tracer, rate, Weight{-1}()),
            routed...,
        )
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
    binding = MortalityParameterBinding(definition, process_id(named))
    return process_fluxes(named, topology, binding)
end
