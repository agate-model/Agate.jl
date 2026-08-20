"""Resolved parameter names for one heterotrophic consumption process."""
struct ConsumptionParameterBinding{T}
    maximum_rate::Symbol
    half_saturation::Symbol
    assimilation::Symbol
    temperature::T
end

function ConsumptionParameterBinding(definition::NormalizedModelDefinition, id::Symbol)
    hasproperty(definition.processes, id) || throw(
        ArgumentError("normalized model has no process :$id"),
    )
    named = getproperty(definition.processes, id)
    process = named.process
    process isa Consumption || throw(ArgumentError("process :$id is not a Consumption process"))
    process.formulation isa HeterotrophicConsumption || throw(
        ArgumentError("unsupported consumption formulation $(typeof(process.formulation))"),
    )
    return ConsumptionParameterBinding(
        parameter_name(definition, _parameter_requirement(named, (), :maximum_rate)),
        parameter_name(definition, _parameter_requirement(named, (), :half_saturation)),
        parameter_name(definition, _parameter_requirement(named, (), :assimilation)),
        _temperature_factor_binding(definition, named),
    )
end

"""Realized consumer-by-resource topology for heterotrophic consumption."""
struct ConsumptionTopology{CT,CI,RT,D}
    consumer_tracers::CT
    consumer_indices::CI
    resource_tracers::RT
    unassimilated_target::D
end

struct ConsumptionResourceLossContribution{F,T} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    consumer_index::Int
    consumer_axis::Int
    resource::Symbol
    resource_axis::Int
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    assimilation_parameter::Symbol
    formulation::F
    temperature_factor::T
end

struct ConsumptionConsumerGainContribution{F,T} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    consumer_index::Int
    consumer_axis::Int
    resource::Symbol
    resource_axis::Int
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    assimilation_parameter::Symbol
    formulation::F
    temperature_factor::T
end

struct ConsumptionUnassimilatedContribution{F,T} <: AbstractProcessContribution
    process::Symbol
    target::Symbol
    consumer_index::Int
    consumer_axis::Int
    resource::Symbol
    resource_axis::Int
    maximum_rate_parameter::Symbol
    half_saturation_parameter::Symbol
    assimilation_parameter::Symbol
    formulation::F
    temperature_factor::T
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
    return ConsumptionTopology(consumer_tracers, consumer_indices, resources, destination)
end

function _consumption_contribution_fields(
    named::NamedProcess,
    topology::ConsumptionTopology,
    binding::ConsumptionParameterBinding,
    consumer_axis::Int,
    resource_axis::Int,
)
    return (
        process_id(named),
        topology.consumer_indices[consumer_axis],
        consumer_axis,
        topology.resource_tracers[resource_axis],
        resource_axis,
        binding.maximum_rate,
        binding.half_saturation,
        binding.assimilation,
        formulation(named.process),
        binding.temperature,
    )
end

function process_contributions(
    named::NamedProcess{P},
    topology::ConsumptionTopology,
    binding::ConsumptionParameterBinding,
) where {P<:Consumption}
    contributions = ()
    for consumer_axis in eachindex(topology.consumer_tracers)
        consumer = topology.consumer_tracers[consumer_axis]
        for resource_axis in eachindex(topology.resource_tracers)
            resource = topology.resource_tracers[resource_axis]
            fields = _consumption_contribution_fields(
                named, topology, binding, consumer_axis, resource_axis
            )
            loss = ConsumptionResourceLossContribution(
                fields[1], resource, fields[2:end]...
            )
            gain = ConsumptionConsumerGainContribution(
                fields[1], consumer, fields[2:end]...
            )
            contributions = (contributions..., loss, gain)
            if !isnothing(topology.unassimilated_target)
                routed = ConsumptionUnassimilatedContribution(
                    fields[1], topology.unassimilated_target, fields[2:end]...
                )
                contributions = (contributions..., routed)
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
) where {P<:Consumption}
    topology = realize_process_topology(named, layout, context)
    binding = ConsumptionParameterBinding(definition, process_id(named))
    return process_contributions(named, topology, binding)
end

struct ConsumptionTerm{F,CI,CA,R,RA,M,K,A,Effect,T}
    formulation::F
    temperature::T
end

@inline _consumption_effect(::Val{:resource_loss}, rate, assimilation) = -rate
@inline _consumption_effect(::Val{:consumer_gain}, rate, assimilation) = assimilation * rate
@inline _consumption_effect(::Val{:unassimilated}, rate, assimilation) =
    (one(assimilation) - assimilation) * rate

@inline function (term::ConsumptionTerm{F,CI,CA,R,RA,M,K,A,Effect,T})(
    bgc, args
) where {F,CI,CA,R,RA,M,K,A,Effect,T}
    consumer = bgc.tracers.plankton(args, CI)
    resource = getproperty(bgc.tracers, R)(args)
    maximum_rate = @inbounds getproperty(bgc.parameters, M)[CA]
    half_saturation = @inbounds getproperty(bgc.parameters, K)[RA]
    assimilation = @inbounds getproperty(bgc.parameters, A)[CA, RA]
    rate = process_rate(term.formulation, resource, consumer, maximum_rate, half_saturation)
    rate = _apply_factor(term.temperature, bgc, args, rate)
    return _consumption_effect(Val(Effect), rate, assimilation)
end

function _consumption_term(contribution, ::Val{Effect}) where {Effect}
    temperature = _lower_factor(contribution.temperature_factor)
    return ConsumptionTerm{
        typeof(contribution.formulation),
        contribution.consumer_index,
        contribution.consumer_axis,
        contribution.resource,
        contribution.resource_axis,
        contribution.maximum_rate_parameter,
        contribution.half_saturation_parameter,
        contribution.assimilation_parameter,
        Effect,
        typeof(temperature),
    }(contribution.formulation, temperature)
end

_lower_contribution(contribution::ConsumptionResourceLossContribution) =
    _consumption_term(contribution, Val(:resource_loss))
_lower_contribution(contribution::ConsumptionConsumerGainContribution) =
    _consumption_term(contribution, Val(:consumer_gain))
_lower_contribution(contribution::ConsumptionUnassimilatedContribution) =
    _consumption_term(contribution, Val(:unassimilated))
