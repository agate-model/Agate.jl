"""Group-level editing helpers for interaction matrices.

Agate accepts several forms for overriding interaction matrices (see
`normalize_interaction_overrides`). For *role-aware* matrices with axes
`(:consumer, :prey)`, a particularly convenient form is a small **group-by-group
block matrix**.

`InteractionBlocks` stores:

- the consumer-group ordering,
- the prey-group ordering,
- the underlying group-block matrix.

It is a host-side convenience object: it is normalized during construction into
the canonical rectangular consumer-by-prey matrix used at runtime.

The editing helpers `set_block!`, `scale_block!`, and `forbid_link!` let users
toggle group links without manually tracking group indices.
"""

"""A group-by-group block matrix for a consumer-by-prey interaction.

Fields
------
- `consumer_groups`: tuple of group symbols defining the **row** ordering
- `prey_groups`: tuple of group symbols defining the **column** ordering
- `B`: the group-block matrix (size `length(consumer_groups) x length(prey_groups)`)
"""
struct InteractionBlocks{CG,PG,M}
    consumer_groups::CG
    prey_groups::PG
    B::M
end

"""Normalize `x` into either `nothing`, a `Tuple{Vararg{Symbol}}`, or return `x`.

This helper is intentionally conservative: it is meant for group-symbol inputs.
"""
@inline function _as_group_tuple(x)
    x === nothing && return nothing
    x isa Symbol && return (x,)
    x isa Tuple && return x
    x isa AbstractVector{Symbol} && return Tuple(x)
    return x
end

"""Convenience constructor for role membership by group symbols.

```julia
roles = roles_from_groups(consumers = :Z, prey = (:P, :Z))
```

This returns the `NamedTuple` accepted by `construct(...; roles=roles)`.
"""
@inline function roles_from_groups(; consumers=nothing, prey=nothing)
    return (consumers=_as_group_tuple(consumers), prey=_as_group_tuple(prey))
end

"""Create a consumer-by-prey group-block interaction object.

The group order is taken from (in order of precedence):

1. `roles` (when provided), otherwise
2. `default_roles(factory)`.

If either role axis is `nothing`, the axis uses `required_groups(factory)`.

The returned `InteractionBlocks` can be edited with `set_block!` and passed
directly as an interaction override.
"""
function interaction_blocks(
    factory::AbstractBGCFactory;
    roles=nothing,
    init=0,
    ::Type{T}=Float64,
) where {T}
    roles_resolved = isnothing(roles) ? default_roles(factory) : roles

    # Roles are expected to be group symbols (or `nothing`).
    consumers = _as_group_tuple(getproperty(roles_resolved, :consumers))
    prey = _as_group_tuple(getproperty(roles_resolved, :prey))

    consumers === nothing && (consumers = required_groups(factory))
    prey === nothing && (prey = required_groups(factory))

    nc = length(consumers)
    np = length(prey)

    B = init == 0 ? zeros(T, nc, np) : fill(T(init), nc, np)
    return InteractionBlocks(consumers, prey, B)
end

@inline function _group_slot(groups, g::Symbol)
    for (i, s) in pairs(groups)
        s === g && return i
    end
    return 0
end

"""Set a consumer-by-prey group block to a value."""
function set_block!(
    blocks::InteractionBlocks;
    consumer_group::Symbol,
    prey_group::Symbol,
    value,
)
    ic = _group_slot(blocks.consumer_groups, consumer_group)
    ic == 0 && throw(ArgumentError("Unknown consumer_group $consumer_group"))
    ip = _group_slot(blocks.prey_groups, prey_group)
    ip == 0 && throw(ArgumentError("Unknown prey_group $prey_group"))
    @inbounds blocks.B[ic, ip] = value
    return blocks
end

"""Scale a consumer-by-prey group block by a factor."""
function scale_block!(
    blocks::InteractionBlocks;
    consumer_group::Symbol,
    prey_group::Symbol,
    factor,
)
    ic = _group_slot(blocks.consumer_groups, consumer_group)
    ic == 0 && throw(ArgumentError("Unknown consumer_group $consumer_group"))
    ip = _group_slot(blocks.prey_groups, prey_group)
    ip == 0 && throw(ArgumentError("Unknown prey_group $prey_group"))
    @inbounds blocks.B[ic, ip] *= factor
    return blocks
end

"""Convenience: set a consumer-prey group link to zero."""
@inline function forbid_link!(
    blocks::InteractionBlocks;
    consumer_group::Symbol,
    prey_group::Symbol,
)
    return set_block!(
        blocks;
        consumer_group=consumer_group,
        prey_group=prey_group,
        value=zero(eltype(blocks.B)),
    )
end
