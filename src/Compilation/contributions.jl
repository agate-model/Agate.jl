"""Setup-time description of one process effect on one concrete tracer."""
abstract type AbstractProcessContribution end

@inline contribution_target(contribution::AbstractProcessContribution) = contribution.target

function _parameter_requirement(
    named::NamedProcess, path::Tuple, slot::Symbol; qualifier::Union{Nothing,NamedTuple}=nothing
)
    matches = filter(parameter_requirements(named)) do requirement
        requirement.identity.path == path &&
        requirement.identity.slot === slot &&
        (isnothing(qualifier) || requirement.identity.qualifier == qualifier)
    end
    length(matches) == 1 || throw(
        ArgumentError(
            isnothing(qualifier) ?
            "process :$(process_id(named)) must declare exactly one $path/$slot parameter requirement" :
            "process :$(process_id(named)) must declare exactly one $path/$slot parameter requirement with qualifier $qualifier",
        ),
    )
    return only(matches)
end

function _scalar_component_target(layout::ComponentLayout, component::Symbol)
    hasproperty(layout.component_tracers, component) || throw(
        ArgumentError("unknown scalar component :$component"),
    )
    tracers = getproperty(layout.component_tracers, component)
    length(tracers) == 1 || throw(
        ArgumentError("component :$component must realize to one tracer"),
    )
    return only(tracers)
end

function _realize_population_classes(
    named::NamedProcess,
    populations::Tuple,
    layout::ComponentLayout,
    context::CommunityContext,
)
    tracer_values = ()
    index_values = ()

    for population in populations
        hasproperty(layout.component_tracers, population) || throw(
            ArgumentError("process :$(named.id) references unrealized population :$population"),
        )
        tracers = getproperty(layout.component_tracers, population)
        indices = Tuple(findfirst(==(tracer), context.plankton_symbols) for tracer in tracers)
        any(isnothing, indices) && throw(
            ArgumentError(
                "process :$(named.id) population :$population realizes tracers absent from the current runtime community",
            ),
        )
        tracer_values = (tracer_values..., tracers...)
        index_values = (index_values..., Int.(indices)...)
    end

    return tracer_values, index_values
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
        Tuple(
            contribution for contribution in contributions if
            contribution_target(contribution) === target
        )
    end
    return NamedTuple{targets}(grouped)
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

"""Derive all typed process contributions for a normalized model."""
function model_contributions(
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext,
)
    contributions = ()
    for named in values(definition.processes)
        contributions = (
            contributions...,
            process_contributions(named, definition, layout, context)...,
        )
    end
    return contributions
end

"""Compile a normalized model into one static equation per requested concrete tracer."""
function compile_model_tendencies(
    definition::NormalizedModelDefinition,
    layout::ComponentLayout,
    context::CommunityContext;
    target_order::Tuple,
)
    contributions = model_contributions(definition, layout, context)
    grouped = group_contributions(contributions; target_order)
    return compile_tendencies(grouped)
end
