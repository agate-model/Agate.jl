"""Resolved parameter names needed to compile one preferential-grazing process."""
struct GrazingParameterBinding
    maximum_rate::Symbol
    half_saturation::Symbol
    palatability::Symbol
    assimilation::Symbol
end

"""Resolve grazing parameter names from normalized semantic requirement bindings."""
function GrazingParameterBinding(definition::NormalizedModelDefinition, id::Symbol)
    hasproperty(definition.processes, id) || throw(
        ArgumentError("normalized model has no process :$id"),
    )
    named = getproperty(definition.processes, id)
    process = named.process
    process isa Grazing || throw(ArgumentError("process :$id is not a Grazing process"))
    process.formulation isa PreferentialGrazing || throw(
        ArgumentError("unsupported grazing formulation $(typeof(process.formulation))"),
    )

    maximum_rate = parameter_name(
        definition, _parameter_requirement(named, (), :maximum_rate)
    )
    half_saturation = parameter_name(
        definition, _parameter_requirement(named, (), :half_saturation)
    )
    palatability = parameter_name(
        definition, _parameter_requirement(named, (), :palatability)
    )
    assimilation = parameter_name(
        definition, _parameter_requirement(named, (), :assimilation)
    )
    return GrazingParameterBinding(
        maximum_rate, half_saturation, palatability, assimilation
    )
end

"""Realized preferential-grazing topology over process-local interaction axes."""
struct GrazingTopology{CT,CI,RT,RI,D}
    consumer_tracers::CT
    consumer_indices::CI
    resource_tracers::RT
    resource_indices::RI
    unassimilated_target::D
end

"""Resource loss from one consumer-by-resource grazing-rate element."""
struct GrazingResourceLossContribution{F} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    consumer::Symbol
    consumer_index::Int
    consumer_axis::Int
    resource_index::Int
    resource_axis::Int
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    palatability_parameter::Symbol
    assimilation_parameter::Symbol
    formulation::F
end

"""Assimilated consumer gain from one consumer-by-resource grazing-rate element."""
struct GrazingConsumerGainContribution{F} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    resource::Symbol
    consumer_index::Int
    consumer_axis::Int
    resource_index::Int
    resource_axis::Int
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    palatability_parameter::Symbol
    assimilation_parameter::Symbol
    formulation::F
end

"""Unassimilated product gain from one consumer-by-resource grazing-rate element."""
struct GrazingUnassimilatedContribution{F} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    consumer::Symbol
    resource::Symbol
    consumer_index::Int
    consumer_axis::Int
    resource_index::Int
    resource_axis::Int
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    palatability_parameter::Symbol
    assimilation_parameter::Symbol
    formulation::F
end

"""Resolve a preferential-grazing process onto process-owned consumer/resource axes."""
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
    return GrazingTopology(
        consumer_tracers,
        consumer_indices,
        resource_tracers,
        resource_indices,
        destination,
    )
end

function _grazing_contribution_fields(
    named::NamedProcess,
    topology::GrazingTopology,
    binding::GrazingParameterBinding,
    consumer_axis::Int,
    resource_axis::Int,
)
    return (
        process_id(named),
        topology.consumer_indices[consumer_axis],
        consumer_axis,
        topology.resource_indices[resource_axis],
        resource_axis,
        binding.maximum_rate,
        binding.half_saturation,
        binding.palatability,
        binding.assimilation,
        formulation(named.process),
    )
end

"""Derive coupled resource, consumer, and unassimilated grazing contributions."""
function process_contributions(
    named::NamedProcess{P}, topology::GrazingTopology, binding::GrazingParameterBinding
) where {P<:Grazing}
    length(topology.consumer_tracers) == length(topology.consumer_indices) || throw(
        ArgumentError("grazing topology consumer tracer and index counts must match"),
    )
    length(topology.resource_tracers) == length(topology.resource_indices) || throw(
        ArgumentError("grazing topology resource tracer and index counts must match"),
    )

    contributions = ()
    for consumer_axis in eachindex(topology.consumer_tracers)
        consumer = topology.consumer_tracers[consumer_axis]
        for resource_axis in eachindex(topology.resource_tracers)
            resource = topology.resource_tracers[resource_axis]
            fields = _grazing_contribution_fields(
                named, topology, binding, consumer_axis, resource_axis
            )
            loss = GrazingResourceLossContribution(
                fields[1], resource, consumer, fields[2:end]...
            )
            gain = GrazingConsumerGainContribution(
                fields[1], consumer, resource, fields[2:end]...
            )
            contributions = (contributions..., loss, gain)

            if !isnothing(topology.unassimilated_target)
                unassimilated = GrazingUnassimilatedContribution(
                    fields[1],
                    topology.unassimilated_target,
                    consumer,
                    resource,
                    fields[2:end]...,
                )
                contributions = (contributions..., unassimilated)
            end
        end
    end
    return contributions
end

function process_contributions(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Grazing}
    topology = realize_process_topology(named, layout, context)
    binding = GrazingParameterBinding(definition, process_id(named))
    return process_contributions(named, topology, binding)
end

struct GrazingTerm{F,CI,CA,RI,RA,M,K,PR,AR,Effect}
    formulation::F
end

function _interaction_runtime_parameter_name(name::Symbol)
    text = String(name)
    suffix = "_matrix"
    return endswith(text, suffix) ? Symbol(text[1:(end - length(suffix))]) : name
end

@inline _grazing_partition_effect(::Val{:consumer_gain}, rate, assimilation) =
    assimilation * rate
@inline _grazing_partition_effect(::Val{:unassimilated}, rate, assimilation) =
    (one(assimilation) - assimilation) * rate

@inline function (term::GrazingTerm{F,CI,CA,RI,RA,M,K,PR,AR,Effect})(
    bgc, args
) where {F,CI,CA,RI,RA,M,K,PR,AR,Effect}
    consumer = bgc.tracers.plankton(args, CI)
    resource = bgc.tracers.plankton(args, RI)
    maximum_rate = @inbounds getproperty(bgc.parameters, M)[CI]
    half_saturation = @inbounds getproperty(bgc.parameters, K)[CI]
    interactions = bgc.parameters.interactions
    palatability = @inbounds getproperty(interactions, PR)[CA, RA]
    rate = process_rate(
        term.formulation,
        resource,
        consumer,
        maximum_rate,
        half_saturation,
        palatability,
    )
    Effect === :resource_loss && return -rate
    assimilation = @inbounds getproperty(interactions, AR)[CA, RA]
    return _grazing_partition_effect(Val(Effect), rate, assimilation)
end

function _grazing_term(contribution, ::Val{Effect}) where {Effect}
    F = typeof(contribution.formulation)
    palatability_runtime = _interaction_runtime_parameter_name(
        contribution.palatability_parameter
    )
    assimilation_runtime = _interaction_runtime_parameter_name(
        contribution.assimilation_parameter
    )
    return GrazingTerm{
        F,
        contribution.consumer_index,
        contribution.consumer_axis,
        contribution.resource_index,
        contribution.resource_axis,
        contribution.maximum_rate_parameter,
        contribution.half_saturation_parameter,
        palatability_runtime,
        assimilation_runtime,
        Effect,
    }(contribution.formulation)
end

_lower_contribution(contribution::GrazingResourceLossContribution) =
    _grazing_term(contribution, Val(:resource_loss))
_lower_contribution(contribution::GrazingConsumerGainContribution) =
    _grazing_term(contribution, Val(:consumer_gain))
_lower_contribution(contribution::GrazingUnassimilatedContribution) =
    _grazing_term(contribution, Val(:unassimilated))
