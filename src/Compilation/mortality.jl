"""Resolved parameter names needed to compile one mortality process instance."""
struct MortalityParameterBinding{R}
    rate::Symbol
    routing_fraction::Union{Nothing,Symbol}
    routing::R
end

MortalityParameterBinding(
    rate::Symbol; routing_fraction::Union{Nothing,Symbol}=nothing
) = MortalityParameterBinding(rate, routing_fraction, nothing)

"""Resolve mortality parameter names from normalized semantic requirement bindings."""
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

"""Loss of one realized population class through one mortality rate element."""
struct MortalityLossContribution{F<:AbstractFormulation} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    source::Symbol
    population_index::Int
    rate_parameter::Symbol
    formulation::F
end

"""Routed gain from one mortality rate element into a product target."""
struct MortalityRoutingContribution{F<:AbstractFormulation} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    source::Symbol
    population_index::Int
    rate_parameter::Symbol
    routing_fraction_parameter::Symbol
    route::Symbol
    formulation::F
end


struct MortalityRoutedProductContribution{F<:AbstractFormulation,Q} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    source::Symbol
    population_index::Int
    rate_parameter::Symbol
    routing_fraction_parameter::Symbol
    ratio_parameter::Q
    route::Symbol
    formulation::F
end

"""Resolve a named mortality process onto the current concrete class layout.

The `ComponentLayout` supplies logical-component tracer identities while the current
`CommunityContext` supplies the flattened plankton indices used by runtime parameter
vectors. The two views are checked for identical realized tracer identities.
"""
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
    if isnothing(process.routing)
        return nothing
    end
    isnothing(binding.routing_fraction) && throw(
        ArgumentError("routed mortality requires a routing-fraction parameter binding"),
    )
    return nothing
end


function _mortality_routed_products(
    named::NamedProcess,
    topology::DOMPOMMortalityTopology,
    binding::MortalityParameterBinding,
    tracer::Symbol,
    population_index::Int,
)
    routing_binding = binding.routing
    isnothing(routing_binding) && throw(
        ArgumentError("DOM/POM-routed mortality requires a routing parameter binding"),
    )
    contributions = ()
    for route in (:DOM, :POM)
        targets = getproperty(topology.routing, route)
        for currency in keys(targets)
            target = getproperty(targets, currency)
            ratio = _routing_ratio_parameter(topology.routing, routing_binding, currency)
            contribution = MortalityRoutedProductContribution(
                process_id(named),
                target,
                tracer,
                population_index,
                binding.rate,
                routing_binding.fraction,
                ratio,
                route,
                formulation(named.process),
            )
            contributions = (contributions..., contribution)
        end
    end
    return contributions
end

"""Derive typed concrete-tracer contributions for one realized mortality process."""
function process_contributions(
    named::NamedProcess{P}, topology::MortalityTopology, binding::MortalityParameterBinding
) where {P<:Mortality}
    length(topology.population_tracers) == length(topology.population_indices) || throw(
        ArgumentError("mortality topology tracer and index counts must match"),
    )
    process = named.process
    _validate_binding(process, binding)
    form = formulation(process)
    routing_fraction = isnothing(process.routing) ? nothing : something(binding.routing_fraction)
    contributions = ()

    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        population_index = topology.population_indices[i]
        loss = MortalityLossContribution(
            process_id(named), tracer, tracer, population_index, binding.rate, form
        )
        contributions = (contributions..., loss)

        if !isnothing(topology.retained_target)
            retained = MortalityRoutingContribution(
                process_id(named),
                topology.retained_target,
                tracer,
                population_index,
                binding.rate,
                routing_fraction,
                :retained,
                form,
            )
            exported = MortalityRoutingContribution(
                process_id(named),
                topology.exported_target,
                tracer,
                population_index,
                binding.rate,
                routing_fraction,
                :exported,
                form,
            )
            contributions = (contributions..., retained, exported)
        end
    end
    return contributions
end


function process_contributions(
    named::NamedProcess{P},
    topology::DOMPOMMortalityTopology,
    binding::MortalityParameterBinding,
) where {P<:Mortality}
    length(topology.population_tracers) == length(topology.population_indices) || throw(
        ArgumentError("mortality topology tracer and index counts must match"),
    )
    contributions = ()
    for i in eachindex(topology.population_tracers)
        tracer = topology.population_tracers[i]
        population_index = topology.population_indices[i]
        loss = MortalityLossContribution(
            process_id(named),
            tracer,
            tracer,
            population_index,
            binding.rate,
            formulation(named.process),
        )
        routed = _mortality_routed_products(
            named, topology, binding, tracer, population_index
        )
        contributions = (contributions..., loss, routed...)
    end
    return contributions
end

function process_contributions(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Mortality}
    topology = realize_process_topology(named, layout, context)
    binding = MortalityParameterBinding(definition, process_id(named))
    return process_contributions(named, topology, binding)
end

struct MortalityLossTerm{F,I,R}
    formulation::F
end

struct MortalityRoutingTerm{F,I,R,Q,Route}
    formulation::F
end


struct MortalityRoutedProductTerm{F,I,R,Q,S,Route}
    formulation::F
end

@inline function (term::MortalityLossTerm{F,I,R})(bgc, args) where {F,I,R}
    biomass = bgc.tracers.plankton(args, I)
    coefficient = @inbounds getproperty(bgc.parameters, R)[I]
    return -process_rate(term.formulation, biomass, coefficient)
end

@inline _routing_weight(::Val{:retained}, fraction) = one(fraction) - fraction
@inline _routing_weight(::Val{:exported}, fraction) = fraction

@inline function (term::MortalityRoutingTerm{F,I,R,Q,Route})(bgc, args) where {F,I,R,Q,Route}
    biomass = bgc.tracers.plankton(args, I)
    coefficient = @inbounds getproperty(bgc.parameters, R)[I]
    fraction = getproperty(bgc.parameters, Q)
    rate = process_rate(term.formulation, biomass, coefficient)
    return _routing_weight(Val(Route), fraction) * rate
end

function _lower_contribution(contribution::MortalityLossContribution{F}) where {F}
    return MortalityLossTerm{
        F,contribution.population_index,contribution.rate_parameter
    }(contribution.formulation)
end

function _lower_contribution(contribution::MortalityRoutingContribution{F}) where {F}
    return MortalityRoutingTerm{
        F,
        contribution.population_index,
        contribution.rate_parameter,
        contribution.routing_fraction_parameter,
        contribution.route,
    }(contribution.formulation)
end


@inline function (term::MortalityRoutedProductTerm{F,I,R,Q,S,Route})(
    bgc, args
) where {F,I,R,Q,S,Route}
    biomass = bgc.tracers.plankton(args, I)
    coefficient = @inbounds getproperty(bgc.parameters, R)[I]
    rate = process_rate(term.formulation, biomass, coefficient)
    fraction = getproperty(bgc.parameters, Q)
    return _organic_route_weight(Val(Route), fraction) *
           _stoichiometric_weight(bgc, Val(S), fraction) * rate
end

function _lower_contribution(contribution::MortalityRoutedProductContribution{F,Q}) where {F,Q}
    return MortalityRoutedProductTerm{
        F,
        contribution.population_index,
        contribution.rate_parameter,
        contribution.routing_fraction_parameter,
        contribution.ratio_parameter,
        contribution.route,
    }(contribution.formulation)
end
