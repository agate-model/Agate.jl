# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Shapes
# ----------------------------------------------------------------------------

const _SHAPES = (:scalar, :vector, :matrix)

@inline function _check_shape(shape::Symbol)
    shape in _SHAPES || throw(ArgumentError("Unknown parameter shape $(shape). Expected one of $(_SHAPES)."))
    return nothing
end

# ----------------------------------------------------------------------------
# Provider values (CPU-only)
# ----------------------------------------------------------------------------

"""Union type for scalar items used inside vector-shaped parameters.

Scalar items may be:
- `Number` / `Bool` literals
- `AbstractParamDef` allometric definitions

`nothing` is intentionally unsupported for vector items (strict-by-default).
"""
const ScalarItem = Union{Number,Bool,AbstractParamDef}

"""Group-level vector provider with a fixed group set.

`GroupVec` stores one scalar item per *group* (e.g. `:P`, `:Z`) in a canonical order.
It is intended for parameters that are constant within each group, independent of
size classes / PFT index.

Guarantees
----------
- The group set is fixed by `groups`.
- The provider is **complete**: every group in `groups` has an explicit value.
- Group values are applied uniformly to all PFTs belonging to that group.

Notes
-----
`GroupVec` is a CPU-side provider stored in the registry. During model construction,
it is expanded to a dense per-PFT vector with minimal branching.
"""
struct GroupVec{N}
    groups::NTuple{N,Symbol}
    items::NTuple{N,ScalarItem}
end

"""\
    GroupVec(groups::NTuple{N,Symbol}, values::NamedTuple) -> GroupVec{N}

Create a `GroupVec` from a complete `NamedTuple` mapping group symbols to scalar items.

The `NamedTuple` must contain **exactly** the keys in `groups` (order is irrelevant).
"""
function GroupVec(groups::NTuple{N,Symbol}, values::NamedTuple) where {N}
    # Strict key validation (typo protection).
    got = Tuple(keys(values))
    missing = Symbol[]
    extra = Symbol[]
    for g in groups
        (g in got) || push!(missing, g)
    end
    for k in got
        (k in groups) || push!(extra, k)
    end
    if !isempty(missing) || !isempty(extra)
        msg = "GroupVec keys must exactly match groups=$(groups)."
        !isempty(missing) && (msg *= " Missing: $(missing).")
        !isempty(extra) && (msg *= " Unknown: $(extra).")
        throw(ArgumentError(msg))
    end

    raw_items = ntuple(i -> begin
        g = groups[i]
        getproperty(values, g)
    end, N)

    # Explicit `nothing` check for actionable errors.
    for i in 1:N
        raw_items[i] === nothing && throw(ArgumentError(
            "GroupVec entries must be explicitly provided for all groups; got `nothing` for group :$(groups[i]).",
        ))
    end

    items = ntuple(i -> (raw_items[i]::ScalarItem), N)
    return GroupVec{N}(groups, items)
end

"""    GroupVec(groups::NTuple{N,Symbol}; kwargs...) -> GroupVec{N}

Keyword form for constructing a complete `GroupVec`.
"""
GroupVec(groups::NTuple{N,Symbol}; kwargs...) where {N} = GroupVec(groups, (; kwargs...))

"""Allowed provider value types stored in a `ParamSpec`.

Provider values are normalized at registry update/extension time.
"""
const ProviderValue = Union{
    Nothing,
    Number,
    Bool,
    AbstractVector,
    AbstractMatrix,
    AbstractParamDef,
    GroupVec,
}

# ----------------------------------------------------------------------------
# Registry types
# ----------------------------------------------------------------------------

"""Single parameter specification.

Fields
------
- `name`: parameter key (Symbol).
- `shape`: one of `:scalar`, `:vector`, `:matrix`.
- `value_kind`: `:real` or `:bool`.
- `doc`: documentation string.
- `provider`: normalized provider value (or `nothing` meaning required).
"""
struct ParamSpec
    name::Symbol
    shape::Symbol
    value_kind::Symbol
    doc::String
    provider::ProviderValue
end

"""Per-model registry (CPU-only)."""
struct ParamRegistry
    specs::Vector{ParamSpec}
    index::Dict{Symbol,Int}
end

"""Create a `ParamRegistry` from a list of `ParamSpec`.

The registry builds a name->index map for fast validation and lookup.
"""
function ParamRegistry(specs::Vector{ParamSpec})
    idx = Dict{Symbol,Int}()
    for (i, s) in pairs(specs)
        haskey(idx, s.name) && throw(ArgumentError("Duplicate ParamSpec for :$(s.name) in registry"))
        idx[s.name] = i
    end
    return ParamRegistry(specs, idx)
end

"""Return the parameter registry for a factory/model."""
function parameter_registry end

@inline function _lookup_index(reg::ParamRegistry, name::Symbol)
    return get(reg.index, name, 0)
end

"""Lookup a `ParamSpec` by name."""
@inline function lookup(reg::ParamRegistry, name::Symbol)
    i = _lookup_index(reg, name)
    i == 0 && return nothing
    return reg.specs[i]
end

@inline function _require_spec(reg::ParamRegistry, name::Symbol)
    i = _lookup_index(reg, name)
    i == 0 && throw(ArgumentError("No ParamSpec found for :$name."))
    return reg.specs[i]
end
