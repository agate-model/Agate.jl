"""ClassRef handles for size-structured plankton communities.

Agate's public API uses **group symbols** (for example `:P`, `:Z` instead of e.g. `:P1`, `:P2`) to describe
roles and to configure interaction matrices. 
Individual size classes within a group are addressed using `ClassRef` handles which include the group symbol and an index of the group which starts at 1.

Fields
------
- `group`: group symbol (e.g. `:P`, `:Z`)
- `i`: index of the class within the group (starting at 1)
"""
struct ClassRef
    group::Symbol
    i::Int
end

"""Construct a `ClassRef`.

A thin wrapper around the `ClassRef` constructor for convenient syntax.

```julia
cref = class(:P, 1)
```
"""
@inline class(group::Symbol, i::Integer) = ClassRef(group, Int(i))

"""Return the number of classes in `group` within a parsed community `community_context`."""
function class_count(community_context, group::Symbol)::Int
    idx = get(community_context.group_indices, group, nothing)
    idx === nothing && throw(ArgumentError("Unknown group symbol $group"))
    return length(idx)
end

"""Resolve `cref` to the global plankton-class index within `community_context`.

The returned index starts at 1 and refers to the flattened plankton ordering
underlying both `community_context.plankton_symbols` and `community_context.group_symbols`.

This is intended for host-side utilities such as initial conditions and diagnostics.
"""
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
