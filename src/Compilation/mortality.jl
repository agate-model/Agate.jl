"""Resolved parameter names needed to compile one mortality process instance."""
struct MortalityParameterBinding
    rate::Symbol
    routing_fraction::Union{Nothing,Symbol}
end

MortalityParameterBinding(
    rate::Symbol; routing_fraction::Union{Nothing,Symbol}=nothing
) = MortalityParameterBinding(rate, routing_fraction)

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
    isnothing(named.process.routing) && return MortalityParameterBinding(rate)
    routing_fraction = parameter_name(
        definition, _parameter_requirement(named, (:routing,), :export_fraction)
    )
    return MortalityParameterBinding(rate; routing_fraction)
end

"""Realized mortality topology over concrete population classes and routing targets."""
struct MortalityTopology{TR,IX,R,E}
    population_tracers::TR
    population_indices::IX
    retained_target::R
    exported_target::E
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
    routing.formulation isa PartitionRouting || throw(
        ArgumentError("unsupported mortality product routing $(typeof(routing.formulation))"),
    )
    retained = _scalar_component_target(layout, routing.retained)
    exported = _scalar_component_target(layout, routing.exported)
    return MortalityTopology(tracer_values, index_values, retained, exported)
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
