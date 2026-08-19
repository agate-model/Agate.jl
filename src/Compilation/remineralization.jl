"""Resolved parameter name needed to compile one linear remineralization process."""
struct RemineralizationParameterBinding
    rate::Symbol
end

"""Resolve a remineralization rate parameter from normalized semantic bindings."""
function RemineralizationParameterBinding(
    definition::NormalizedModelDefinition, id::Symbol
)
    hasproperty(definition.processes, id) || throw(
        ArgumentError("normalized model has no process :$id"),
    )
    named = getproperty(definition.processes, id)
    process = named.process
    process isa Remineralization || throw(
        ArgumentError("process :$id is not a Remineralization process"),
    )
    process.formulation isa LinearRemineralization || throw(
        ArgumentError("unsupported remineralization formulation $(typeof(process.formulation))"),
    )
    rate = parameter_name(definition, _parameter_requirement(named, (), :rate))
    return RemineralizationParameterBinding(rate)
end

"""Realized linear remineralization topology over scalar source and destination tracers."""
struct RemineralizationTopology{S,D}
    source_tracers::S
    destination_target::D
end

"""Loss from one realized remineralization source."""
struct RemineralizationSourceLossContribution{F} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    source::Symbol
    rate_parameter::Symbol
    formulation::F
end

"""Destination gain coupled to one realized remineralization source rate."""
struct RemineralizationDestinationGainContribution{F} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    source::Symbol
    rate_parameter::Symbol
    formulation::F
end

"""Resolve linear remineralization onto scalar source tracers and one destination."""
function realize_process_topology(
    named::NamedProcess{P}, layout::ComponentLayout, context::CommunityContext
) where {P<:Remineralization}
    process = named.process
    process.formulation isa LinearRemineralization || throw(
        ArgumentError("unsupported remineralization formulation $(typeof(process.formulation))"),
    )
    length(process.destinations) == 1 || throw(
        ArgumentError("linear remineralization currently requires one destination"),
    )
    sources = Tuple(_scalar_component_target(layout, source) for source in process.sources)
    destination = _scalar_component_target(layout, only(process.destinations))
    return RemineralizationTopology(sources, destination)
end

"""Derive coupled source-loss and destination-gain remineralization contributions."""
function process_contributions(
    named::NamedProcess{P},
    topology::RemineralizationTopology,
    binding::RemineralizationParameterBinding,
) where {P<:Remineralization}
    process = named.process
    length(topology.source_tracers) == length(process.sources) || throw(
        ArgumentError("remineralization topology source count must match process sources"),
    )
    contributions = ()
    for source in topology.source_tracers
        fields = (process_id(named), source, binding.rate, process.formulation)
        loss = RemineralizationSourceLossContribution(fields[1], source, fields[2:end]...)
        gain = RemineralizationDestinationGainContribution(
            fields[1], topology.destination_target, fields[2:end]...
        )
        contributions = (contributions..., loss, gain)
    end
    return contributions
end

function process_contributions(
    named::NamedProcess{P},
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
) where {P<:Remineralization}
    topology = realize_process_topology(named, layout, context)
    binding = RemineralizationParameterBinding(definition, process_id(named))
    return process_contributions(named, topology, binding)
end

struct RemineralizationTerm{F,S,R,Effect}
    formulation::F
end

@inline _remineralization_effect(::Val{:source_loss}, rate) = -rate
@inline _remineralization_effect(::Val{:destination_gain}, rate) = rate

@inline function (term::RemineralizationTerm{F,S,R,Effect})(bgc, args) where {F,S,R,Effect}
    source = getproperty(bgc.tracers, S)(args)
    coefficient = getproperty(bgc.parameters, R)
    rate = process_rate(term.formulation, source, coefficient)
    return _remineralization_effect(Val(Effect), rate)
end

function _remineralization_term(contribution, ::Val{Effect}) where {Effect}
    F = typeof(contribution.formulation)
    return RemineralizationTerm{
        F, contribution.source, contribution.rate_parameter, Effect
    }(contribution.formulation)
end

_lower_contribution(contribution::RemineralizationSourceLossContribution) =
    _remineralization_term(contribution, Val(:source_loss))
_lower_contribution(contribution::RemineralizationDestinationGainContribution) =
    _remineralization_term(contribution, Val(:destination_gain))
