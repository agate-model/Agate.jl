"""Reference one concrete class within a named structured group or component.

`group` may identify a parsed community subgroup or a logical component in a
`ComponentLayout`. `i` is the one-based class ordinal within that identity.
"""
struct ClassRef
    group::Symbol
    i::Int
end

"""Construct a `ClassRef` for host-side class lookup."""
@inline class(group::Symbol, i::Integer) = ClassRef(group, Int(i))

"""Return the number of classes in a parsed community subgroup."""
function class_count(community_context, group::Symbol)::Int
    idx = get(community_context.group_indices, group, nothing)
    idx === nothing && throw(ArgumentError("Unknown group symbol $group"))
    return length(idx)
end

"""Return the number of concrete classes realized by a logical component."""
class_count(layout::ComponentLayout, component::Symbol)::Int =
    component_class_count(layout, component)

"""Resolve a community `ClassRef` to its flattened community-class index."""
function resolve_class(community_context, cref::ClassRef)::Int
    idx = get(community_context.group_indices, cref.group, nothing)
    idx === nothing && throw(ArgumentError("Unknown group symbol $(cref.group)"))
    1 <= cref.i <= length(idx) || throw(
        ArgumentError(
            "Class ordinal $(cref.i) is out of bounds for group $(cref.group) (valid 1:$(length(idx))).",
        ),
    )
    return idx[cref.i]
end

"""Resolve a component `ClassRef` to one physical tracer position.

Pool classes and single-state population classes are unambiguous. Multi-state
populations require an explicit state reference and therefore cannot be resolved
by `ClassRef` alone.
"""
function resolve_class(layout::ComponentLayout, cref::ClassRef)::Int
    nclasses = component_class_count(layout, cref.group)
    1 <= cref.i <= nclasses || throw(
        ArgumentError(
            "Class ordinal $(cref.i) is out of bounds for component $(cref.group) (valid 1:$nclasses).",
        ),
    )
    state_mapping = component_state_tracers(layout, cref.group)
    if state_mapping === nothing
        return component_indices(layout, cref.group)[cref.i]
    end
    length(state_mapping) == 1 || throw(
        ArgumentError(
            "ClassRef for multi-state component $(cref.group) is ambiguous; select a prognostic state explicitly",
        ),
    )
    state = only(keys(state_mapping))
    return state_indices(layout, cref.group, state)[cref.i]
end
