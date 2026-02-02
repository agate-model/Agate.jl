"""Class identity handles for size-structured plankton communities.

Agate's public API uses **group symbols** (for example `:P`, `:Z`) to describe
roles and to configure interaction matrices. Individual size classes within a
group are addressed using `ClassRef` handles rather than "symbol-with-suffix"
names like `:P1`.

`ClassRef` is intentionally tiny and GPU-friendly (isbits). It is meant to be
resolved on the host to concrete indices.
"""

"""A stable reference to a plankton size class.

Fields
------
- `group`: group symbol (e.g. `:P`, `:Z`)
- `i`: 1-based ordinal within the group
"""
struct ClassRef
    group::Symbol
    i::Int
end

"""Construct a `ClassRef`.

```julia
cref = class(:P, 1)
```
"""
@inline class(group::Symbol, i::Integer) = ClassRef(group, Int(i))

"""Return the number of classes in `group` within a parsed community `ctx`."""
function class_count(ctx, group::Symbol)::Int
    idx = get(ctx.group_indices, group, nothing)
    idx === nothing && throw(ArgumentError("Unknown group symbol $group"))
    return length(idx)
end

"""Resolve `cref` to the *global* plankton-class index within `ctx`.

The returned index is 1-based and refers to the flattened plankton ordering
(`ctx.plankton_symbols` / `ctx.group_symbols`).

This is intended for host-side utilities (initial conditions, diagnostics).
"""
function resolve_class(ctx, cref::ClassRef)::Int
    idx = get(ctx.group_indices, cref.group, nothing)
    idx === nothing && throw(ArgumentError("Unknown group symbol $(cref.group)"))
    1 <= cref.i <= length(idx) || throw(
        ArgumentError(
            "Class ordinal $(cref.i) is out of bounds for group $(cref.group) (valid 1:$(length(idx))).",
        ),
    )
    return idx[cref.i]
end
