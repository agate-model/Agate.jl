# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Provider normalization (strict-by-default)
# ----------------------------------------------------------------------------

const _VALUE_KINDS = (:real, :bool)

@inline function _check_value_kind(value_kind::Symbol)
    value_kind in _VALUE_KINDS ||
        throw(ArgumentError("Unknown value_kind $(value_kind). Expected one of $(_VALUE_KINDS)."))
    return nothing
end

@inline _bad_provider(msg::AbstractString) = throw(ArgumentError(msg))

@inline function _normalize_scalar_item(x, value_kind::Symbol)
    x === nothing && return nothing
    _check_value_kind(value_kind)

    if value_kind === :bool
        x isa Bool && return x
        return _bad_provider("Invalid scalar provider item $(typeof(x)). Expected Bool.")
    end

    # :real
    if x isa Number || x isa AbstractParamDef
        return x
    end

    return _bad_provider("Invalid scalar provider item $(typeof(x)). Expected a number or AbstractParamDef.")
end

@inline function _normalize_required_scalar_item(x, value_kind::Symbol)
    x === nothing && return _bad_provider("Vector entries must be explicitly provided; `nothing` is not allowed.")
    return _normalize_scalar_item(x, value_kind)
end

@inline function _check_container_kind(kind::Symbol, elty, shape::Symbol)
    # Keep checks strict but cheap and informative.
    if kind === :bool
        elty <: Bool || _bad_provider("Invalid $(shape) provider element type $(elty). Expected Bool.")
    else
        elty <: Number || _bad_provider("Invalid $(shape) provider element type $(elty). Expected a subtype of Number.")
    end
    return nothing
end

_normalize_vector_provider(::Nothing, value_kind::Symbol) = nothing

function _normalize_vector_provider(x::AbstractVector, value_kind::Symbol)
    _check_value_kind(value_kind)
    _check_container_kind(value_kind, eltype(x), :vector)
    return x
end

function _normalize_vector_provider(x::NamedTuple, ::Symbol)
    return _bad_provider(
        "Vector parameters do not accept per-group NamedTuple maps. " *
        "Use GroupVec(groups; ...) for group-level vectors, " *
        "or provide a full-length AbstractVector / scalar broadcast."
    )
end

function _normalize_vector_provider(x::GroupVec{N}, value_kind::Symbol) where {N}
    _check_value_kind(value_kind)
    items = ntuple(i -> _normalize_required_scalar_item(x.items[i], value_kind), N)
    return GroupVec{N}(x.groups, items)
end

_normalize_vector_provider(::Dict, ::Symbol) = _bad_provider(
    "Vector overrides do not accept Dict inputs (typo-prone). " *
    "Use GroupVec(groups; ...) for full replacement."
)

_normalize_vector_provider(::Function, ::Symbol) = _bad_provider(
    "Vector parameters do not accept function providers. Provide a scalar broadcast, full-length AbstractVector, or GroupVec(groups; ...)."
)

function _normalize_vector_provider(x, value_kind::Symbol)
    # Shape-driven rule: scalar/Bool/allometric inputs for vector parameters broadcast across all PFTs.
    return _normalize_required_scalar_item(x, value_kind)
end

_normalize_matrix_provider(::Nothing, ::Symbol) = nothing

function _normalize_matrix_provider(x::AbstractMatrix, value_kind::Symbol)
    _check_value_kind(value_kind)
    _check_container_kind(value_kind, eltype(x), :matrix)
    return x
end

_normalize_matrix_provider(::Function, ::Symbol) = _bad_provider(
    "Matrix parameters do not accept function providers. Provide an AbstractMatrix, or pass `nothing` to use model defaults.",
)

_normalize_matrix_provider(x, ::Symbol) =
    _bad_provider("Invalid matrix provider $(typeof(x)). Expected AbstractMatrix or `nothing`.")

@inline function _normalize_scalar_provider(x, value_kind::Symbol)
    x === nothing && return nothing
    _check_value_kind(value_kind)

    if x isa Function
        return _bad_provider("Scalar parameters do not accept function providers. Provide a number/Bool.")
    end

    if value_kind === :bool
        x isa Bool && return x
        return _bad_provider("Invalid scalar provider $(typeof(x)). Expected Bool.")
    end

    # :real
    if x isa Number
        return x
    elseif x isa AbstractParamDef
        return _bad_provider(
            "Allometric providers are not valid scalar defaults; use them inside vector-shaped parameters.",
        )
    end

    return _bad_provider("Invalid scalar provider $(typeof(x)). Expected a number.")
end

"""Normalize a user input for a parameter of the declared `shape`.

Normalization happens at registry update/extension time, so runtime resolution can rely on
canonical provider values.
"""
function normalize_provider(shape::Symbol, x, value_kind::Symbol)
    _check_shape(shape)
    _check_value_kind(value_kind)

    if shape === :scalar
        return _normalize_scalar_provider(x, value_kind)
    elseif shape === :vector
        return _normalize_vector_provider(x, value_kind)
    else
        return _normalize_matrix_provider(x, value_kind)
    end
end

# ----------------------------------------------------------------------------
# Public constructors
# ----------------------------------------------------------------------------

"""\
    ParamSpec(name, shape, doc, default; value_kind=:real)

Create a `ParamSpec` and normalize `default` into a canonical provider value.
"""
function ParamSpec(
    name::Symbol,
    shape::Symbol,
    doc::AbstractString,
    default;
    value_kind::Symbol=:real,
)
    _check_shape(shape)
    _check_value_kind(value_kind)
    prov = normalize_provider(shape, default, value_kind)
    return ParamSpec(name, shape, value_kind, String(doc), prov)
end

"""\
    scalar_param(name, doc, default; value_kind=:real) -> ParamSpec

Create a scalar parameter specification.
"""
scalar_param(name::Symbol, doc::AbstractString, default; value_kind::Symbol=:real) =
    ParamSpec(name, :scalar, doc, default; value_kind=value_kind)

"""\
    vector_param(name, doc, default; value_kind=:real) -> ParamSpec

Create a vector parameter specification.
"""
vector_param(name::Symbol, doc::AbstractString, default; value_kind::Symbol=:real) =
    ParamSpec(name, :vector, doc, default; value_kind=value_kind)

"""\
    matrix_param(name, doc, default; value_kind=:real) -> ParamSpec

Create a matrix parameter specification.
"""
matrix_param(name::Symbol, doc::AbstractString, default; value_kind::Symbol=:real) =
    ParamSpec(name, :matrix, doc, default; value_kind=value_kind)
