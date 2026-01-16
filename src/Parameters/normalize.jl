# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Provider normalization (strict-by-default)
# ----------------------------------------------------------------------------

@inline function _normalize_scalar_item(x)
    x === nothing && return nothing

    if x isa Number || x isa Bool || x isa AbstractParamDef
        return x
    end

    throw(ArgumentError("Invalid scalar provider item $(typeof(x)). Expected a number/Bool or AbstractParamDef."))
end

function _normalize_vector_provider(x)
    x === nothing && return nothing

    if x isa AbstractVector
        return x
    elseif x isa NamedTuple
        ks = collect(keys(x))
        vs = values(x)
        items = Vector{ScalarItem}(undef, length(ks))
        for i in eachindex(ks)
            items[i] = _normalize_scalar_item(vs[i])
        end
        return VectorGroupMap(ks, items)
    elseif x isa Dict
        throw(ArgumentError(
            "Vector group maps must be provided as a NamedTuple like (Z=..., P=...). Dict inputs are intentionally unsupported.",
        ))
    elseif x isa Function
        throw(ArgumentError(
            "Vector parameters do not accept function providers. Provide a full-length vector or a per-group NamedTuple like (Z=..., P=...).",
        ))
    else
        # Shape-driven rule: scalar/Bool/allometric inputs for vector parameters broadcast across all PFTs.
        return _normalize_scalar_item(x)
    end
end

function _normalize_matrix_provider(x)
    x === nothing && return nothing

    if x isa AbstractMatrix
        return x
    elseif x isa MatrixFn
        return x
    elseif x isa Function
        throw(ArgumentError(
            "Matrix parameters do not accept bare function providers. Wrap the function in MatrixFn(f; deps=(...)).",
        ))
    else
        throw(ArgumentError("Invalid matrix provider $(typeof(x)). Expected AbstractMatrix or MatrixFn."))
    end
end

"""Normalize a user input for a parameter of the declared `shape`.

Normalization happens at registry update/extension time, so runtime resolution can rely on
the canonical provider values.
"""
function normalize_provider(shape::Symbol, x)
    _check_shape(shape)

    if shape === :scalar
        x === nothing && return nothing

        if x isa Function
            throw(ArgumentError("Scalar parameters do not accept function providers. Provide a number/Bool."))
        elseif x isa Number || x isa Bool
            return x
        elseif x isa AbstractParamDef
            throw(ArgumentError("Allometric providers are not valid scalar defaults; use them inside vector-shaped parameters."))
        end

        throw(ArgumentError("Invalid scalar provider $(typeof(x)). Expected number/Bool."))

    elseif shape === :vector
        return _normalize_vector_provider(x)

    else # :matrix
        return _normalize_matrix_provider(x)
    end
end

# ----------------------------------------------------------------------------
# Public constructors
# ----------------------------------------------------------------------------

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
"""
scalar_param(name::Symbol, doc::AbstractString, default; missing_policy::Symbol=:fail, value_kind::Symbol=:real) =
    ParamSpec(name, :scalar, doc, default; missing_policy=missing_policy, value_kind=value_kind)

"""\
    vector_param(name, doc, default; missing_policy=:fail, value_kind=:real) -> ParamSpec

Create a vector parameter specification.
"""
vector_param(name::Symbol, doc::AbstractString, default; missing_policy::Symbol=:fail, value_kind::Symbol=:real) =
    ParamSpec(name, :vector, doc, default; missing_policy=missing_policy, value_kind=value_kind)

"""\
    matrix_param(name, doc, default; missing_policy=:fail, value_kind=:real) -> ParamSpec

Create a matrix parameter specification.
"""
matrix_param(name::Symbol, doc::AbstractString, default; missing_policy::Symbol=:fail, value_kind::Symbol=:real) =
    ParamSpec(name, :matrix, doc, default; missing_policy=missing_policy, value_kind=value_kind)
