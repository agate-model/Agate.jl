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
- `nothing` (interpreted according to `missing_policy`)
"""
const ScalarItem = Union{Nothing,Number,Bool,AbstractParamDef}

"""Vector provider defined by per-group scalar defaults.

The group map stores scalar items keyed by group symbol (e.g. `:Z`, `:P`).
Per-PFT overrides stored in the community specification take precedence.

Accepted construction inputs:
- `NamedTuple`, e.g. `(Z=1.0, P=0.0)`

`Dict` inputs are intentionally unsupported to keep configuration strict and typo-resistant.
"""
struct VectorGroupMap
    keys::Vector{Symbol}
    items::Vector{ScalarItem}
end

"""    MatrixFn(f; deps=Symbol[])

Derived matrix provider with explicit dependencies.

`f` is called during parameter resolution as `f(ctx, depvals)` where `depvals` is an
`NTuple` of dependency values in the same order as `deps=...`.

`MatrixFn` providers are evaluated during `construct` and must return an `AbstractMatrix`
of size `(n_total, n_total)`.
"""
struct MatrixFn{F,N}
    f::F
    deps::NTuple{N,Symbol}
end

"""    MatrixFn(f; deps=Symbol[])

Create a derived matrix provider with explicit dependencies.

Dependencies are stored as an `NTuple{N,Symbol}` and passed positionally to
`f(ctx, depvals)`.
"""
function MatrixFn(f; deps=Symbol[])
    deps_tuple = Tuple(Symbol(d) for d in deps)
    return MatrixFn{typeof(f), length(deps_tuple)}(f, deps_tuple)
end

"""Return declared dependencies for a provider."""
deps(::Any) = ()
deps(p::MatrixFn) = p.deps

# ----------------------------------------------------------------------------
# Registry types
# ----------------------------------------------------------------------------

"""Single parameter specification.

Fields
------
- `name`: parameter key (Symbol).
- `shape`: one of `:scalar`, `:vector`, `:matrix`.
- `missing_policy`:
  - `:fail`        -> throw on missing/`nothing`
  - `:zero_warn`   -> replace with 0/false and warn (once per parameter)
  - `:zero_silent` -> replace with 0/false silently
- `value_kind`: `:real` or `:bool`.
- `doc`: documentation string.
- `provider`: normalized provider value (or `nothing` meaning required).
"""
struct ParamSpec
    name::Symbol
    shape::Symbol
    missing_policy::Symbol
    value_kind::Symbol
    doc::String
    provider::Any
end

"""\
    ParamSpec(name, shape, doc, default; missing_policy=:fail, value_kind=:real)

Create a `ParamSpec` and normalize `default` into a canonical provider value.
"""
function ParamSpec(
    name::Symbol,
    shape::Symbol,
    doc::AbstractString,
    default;
    missing_policy::Symbol=:fail,
    value_kind::Symbol=:real,
)
    _check_shape(shape)
    prov = default === nothing ? nothing : normalize_provider(shape, default)
    return ParamSpec(name, shape, missing_policy, value_kind, String(doc), prov)
end

"""\
    scalar_param(name, doc, default; missing_policy=:fail, value_kind=:real) -> ParamSpec

Create a scalar parameter specification.

This is the preferred user-facing constructor for new scalar parameters. The internal
registry stores an explicit shape, but callers should not need to pass `:scalar`.
"""
scalar_param(name::Symbol, doc::AbstractString, default; missing_policy::Symbol=:fail, value_kind::Symbol=:real) =
    ParamSpec(name, :scalar, doc, default; missing_policy=missing_policy, value_kind=value_kind)

"""\
    vector_param(name, doc, default; missing_policy=:fail, value_kind=:real) -> ParamSpec

Create a vector parameter specification.

`default` may be a scalar/Bool/allometric definition (broadcast to all PFTs),
a full-length vector, or a per-group `NamedTuple` map (e.g. `(Z=..., P=...)`).
"""
vector_param(name::Symbol, doc::AbstractString, default; missing_policy::Symbol=:fail, value_kind::Symbol=:real) =
    ParamSpec(name, :vector, doc, default; missing_policy=missing_policy, value_kind=value_kind)

"""\
    matrix_param(name, doc, default; missing_policy=:fail, value_kind=:real) -> ParamSpec

Create a matrix parameter specification.

`default` may be a concrete matrix or a `MatrixFn(f; deps=[...])`.
"""
matrix_param(name::Symbol, doc::AbstractString, default; missing_policy::Symbol=:fail, value_kind::Symbol=:real) =
    ParamSpec(name, :matrix, doc, default; missing_policy=missing_policy, value_kind=value_kind)

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
