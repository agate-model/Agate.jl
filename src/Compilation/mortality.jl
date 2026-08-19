"""Setup-time description of one process effect on one concrete tracer."""
abstract type AbstractProcessContribution end

"""Resolved parameter names needed to compile one mortality process instance."""
struct MortalityParameterBinding
    rate::Symbol
    routing_fraction::Union{Nothing,Symbol}
end

MortalityParameterBinding(
    rate::Symbol; routing_fraction::Union{Nothing,Symbol}=nothing
) = MortalityParameterBinding(rate, routing_fraction)

function _mortality_requirement(named::NamedProcess, path::Tuple, slot::Symbol)
    matches = filter(
        requirement ->
            requirement.identity.path == path && requirement.identity.slot === slot,
        parameter_requirements(named),
    )
    length(matches) == 1 || throw(
        ArgumentError(
            "process :$(process_id(named)) must declare exactly one $path/$slot parameter requirement",
        ),
    )
    return only(matches)
end

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

    rate = parameter_name(definition, _mortality_requirement(named, (), :rate))
    isnothing(named.process.routing) && return MortalityParameterBinding(rate)
    routing_fraction = parameter_name(
        definition, _mortality_requirement(named, (:routing,), :export_fraction)
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

@inline contribution_target(contribution::AbstractProcessContribution) = contribution.target

function _scalar_component_target(layout::ComponentLayout, component::Symbol)
    hasproperty(layout.component_tracers, component) || throw(
        ArgumentError("unknown routing component :$component"),
    )
    tracers = getproperty(layout.component_tracers, component)
    length(tracers) == 1 || throw(
        ArgumentError("routing component :$component must realize to one tracer"),
    )
    return only(tracers)
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
    tracer_values = ()
    index_values = ()

    for population in process.populations
        hasproperty(layout.component_tracers, population) || throw(
            ArgumentError("process :$(named.id) references unrealized population :$population"),
        )
        tracers = getproperty(layout.component_tracers, population)
        indices = get(context.group_indices, population, nothing)
        isnothing(indices) && throw(
            ArgumentError(
                "process :$(named.id) population :$population is absent from the current runtime community",
            ),
        )
        runtime_tracers = Tuple(context.plankton_symbols[indices])
        tracers == runtime_tracers || throw(
            ArgumentError(
                "process :$(named.id) population :$population realizes as $tracers but the current runtime community realizes it as $runtime_tracers",
            ),
        )
        tracer_values = (tracer_values..., tracers...)
        index_values = (index_values..., Tuple(indices)...)
    end

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

function _target_order(contributions::Tuple)
    targets = Symbol[]
    for contribution in contributions
        target = contribution_target(contribution)
        target in targets || push!(targets, target)
    end
    return Tuple(targets)
end

function _validate_target_order(contributions::Tuple, target_order::Tuple)
    length(unique(target_order)) == length(target_order) || throw(
        ArgumentError("target_order contains duplicate tracer identities"),
    )
    actual = Set(_target_order(contributions))
    requested = Set(target_order)
    actual == requested || throw(
        ArgumentError("target_order must contain exactly the contribution targets"),
    )
    return nothing
end

"""Group a flat contribution tuple by concrete target tracer.

An explicit `target_order` can be supplied to match an established runtime tracer order.
"""
function group_contributions(contributions::Tuple; target_order=nothing)
    isempty(contributions) && throw(ArgumentError("cannot group an empty contribution tuple"))
    targets = isnothing(target_order) ? _target_order(contributions) : Tuple(target_order)
    _validate_target_order(contributions, targets)

    grouped = ntuple(length(targets)) do i
        target = targets[i]
        Tuple(contribution for contribution in contributions if contribution_target(contribution) === target)
    end
    return NamedTuple{targets}(grouped)
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

struct StaticContributionEquation{T}
    terms::T
end

@inline function _sum_contributions(terms::Tuple{T}, bgc, args) where {T}
    return first(terms)(bgc, args)
end

@inline function _sum_contributions(
    terms::Tuple{T,S,Vararg{Any,N}}, bgc, args
) where {T,S,N}
    return first(terms)(bgc, args) + _sum_contributions(Base.tail(terms), bgc, args)
end

@inline function (equation::StaticContributionEquation)(bgc, x, y, z, t, args...)
    return _sum_contributions(equation.terms, bgc, args)
end

"""Statically lower one target's contribution tuple to a `CompiledEquation`."""
function compile_tendency(contributions::Tuple)
    isempty(contributions) && throw(ArgumentError("cannot compile an empty contribution tuple"))
    target = contribution_target(first(contributions))
    all(contribution -> contribution_target(contribution) === target, contributions) || throw(
        ArgumentError("compile_tendency requires contributions for one target tracer"),
    )
    terms = map(_lower_contribution, contributions)
    return CompiledEquation(StaticContributionEquation(terms))
end

"""Compile every grouped target contribution into a concrete tracer equation."""
compile_tendencies(grouped::NamedTuple) = map(compile_tendency, grouped)
