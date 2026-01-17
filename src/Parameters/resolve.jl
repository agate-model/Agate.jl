# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Resolution helpers
# ----------------------------------------------------------------------------

@inline _is_bool(spec::ParamSpec) = spec.value_kind === :bool

"""Convert a user value to the runtime storage type for a given `value_kind`."""
@inline function coerce_value(::Type{FT}, value_kind::Symbol, x) where {FT<:AbstractFloat}
    if value_kind === :bool
        x isa Bool || throw(ArgumentError("Expected Bool, got $(typeof(x))."))
        return x
    end
    x isa Bool && throw(ArgumentError("Expected a number, got Bool."))
    return FT(x)
end

function _resolve_scalar(::Type{FT}, spec::ParamSpec, ctx) where {FT<:AbstractFloat}
    spec.shape === :scalar || throw(ArgumentError("Internal error: :$(spec.name) is not a scalar spec."))
    p = spec.provider
    isnothing(p) && throw(ArgumentError("Missing required scalar parameter :$(spec.name)."))

    if p isa Number || p isa Bool
        return coerce_value(FT, spec.value_kind, p)
    end

    throw(ArgumentError("Internal error: unsupported scalar provider $(typeof(p)) for :$(spec.name)."))
end

@inline function _resolve_scalar_item(::Type{FT}, spec::ParamSpec, ctx, idx::Int, p) where {FT<:AbstractFloat}
    p === nothing && throw(ArgumentError("Internal error: unexpected `nothing` item for vector :$(spec.name)."))

    if p isa Number || p isa Bool
        return coerce_value(FT, spec.value_kind, p)
    elseif p isa AbstractParamDef
        val = resolve_param(FT, p, ctx.diameters[idx])
        return coerce_value(FT, spec.value_kind, val)
    end

    throw(ArgumentError("Internal error: unsupported vector item provider $(typeof(p)) for :$(spec.name)."))
end

"""Resolve a vector parameter over the full plankton community.

Guarantees
----------
- Vector parameters are strict-by-default: `nothing` is an error.
- Group-level providers (`GroupVec`) are expanded once during construction.

Notes
-----
Per-PFT (index-level) parameter overrides are intentionally unsupported.
"""
function _resolve_vector(::Type{FT}, spec::ParamSpec, ctx) where {FT<:AbstractFloat}
    spec.shape === :vector || throw(ArgumentError("Internal error: :$(spec.name) is not a vector spec."))

    n = ctx.n_total
    is_bool = _is_bool(spec)

    p = spec.provider
    isnothing(p) && throw(ArgumentError("Missing required vector parameter :$(spec.name)."))

    out = is_bool ? Vector{Bool}(undef, n) : Vector{FT}(undef, n)

    if p isa GroupVec
        groups = p.groups
        N = length(groups)
        @inbounds for i in 1:n
            gsym = ctx.group_symbols[i]
            gi = 0
            for j in 1:N
                if groups[j] === gsym
                    gi = j
                    break
                end
            end
            gi == 0 && throw(ArgumentError(
                "Vector :$(spec.name) is defined over groups=$(groups) but the community contains group :$(gsym).",
            ))
            out[i] = _resolve_scalar_item(FT, spec, ctx, i, p.items[gi])
        end

    elseif p isa AbstractVector
        length(p) == n || throw(ArgumentError("Default for :$(spec.name) must have length $n."))
        @inbounds for i in 1:n
            out[i] = _resolve_scalar_item(FT, spec, ctx, i, p[i])
        end

    else
        # Scalar/Bool/allometric broadcasts.
        @inbounds for i in 1:n
            out[i] = _resolve_scalar_item(FT, spec, ctx, i, p)
        end
    end

    return out
end

# ----------------------------------------------------------------------------
# Built-in derived matrices
# ----------------------------------------------------------------------------

using ..Library.Allometry:
    palatability_matrix_allometric,
    allometric_palatability_unimodal_protection,
    assimilation_efficiency_matrix_binary

const _DERIVED_MATRIX_DEPS = Dict{Symbol,NTuple{4,Symbol}}(
    :palatability_matrix => (:can_eat, :optimum_predator_prey_ratio, :specificity, :protection),
)

const _DERIVED_MATRIX_DEPS3 = Dict{Symbol,NTuple{3,Symbol}}(
    :assimilation_matrix => (:can_eat, :can_be_eaten, :assimilation_efficiency),
)

@inline _is_derived_matrix(k::Symbol) = haskey(_DERIVED_MATRIX_DEPS, k) || haskey(_DERIVED_MATRIX_DEPS3, k)

function _derive_matrix(::Type{FT}, name::Symbol, ctx, resolved::Dict{Symbol,Any}) where {FT<:AbstractFloat}
    n = ctx.n_total

    if name === :palatability_matrix
        can_eat = resolved[:can_eat]
        opr = resolved[:optimum_predator_prey_ratio]
        spec = resolved[:specificity]
        prot = resolved[:protection]

        M = palatability_matrix_allometric(FT, ctx.diameters;
            can_eat=can_eat,
            optimum_predator_prey_ratio=opr,
            specificity=spec,
            protection=prot,
            palatability_fn=allometric_palatability_unimodal_protection,
        )
        (size(M, 1) == n && size(M, 2) == n) || throw(ArgumentError(
            "Derived :palatability_matrix must be size ($n,$n).",
        ))
        return M

    elseif name === :assimilation_matrix
        can_eat = resolved[:can_eat]
        can_be_eaten = resolved[:can_be_eaten]
        eff = resolved[:assimilation_efficiency]

        M = assimilation_efficiency_matrix_binary(FT;
            can_eat=can_eat,
            can_be_eaten=can_be_eaten,
            assimilation_efficiency=eff,
        )
        (size(M, 1) == n && size(M, 2) == n) || throw(ArgumentError(
            "Derived :assimilation_matrix must be size ($n,$n).",
        ))
        return M
    end

    throw(ArgumentError("No derived matrix builder registered for :$name."))
end

function _resolve_matrix(::Type{FT}, spec::ParamSpec, ctx, resolved::Dict{Symbol,Any}) where {FT<:AbstractFloat}
    spec.shape === :matrix || throw(ArgumentError("Internal error: :$(spec.name) is not a matrix spec."))

    n = ctx.n_total
    p = spec.provider

    if p isa AbstractMatrix
        (size(p, 1) == n && size(p, 2) == n) || throw(ArgumentError(
            "Matrix :$(spec.name) must be size ($n,$n), got $(size(p)).",
        ))
        # Coerce element type eagerly to keep runtime storage uniform.
        return FT.(p)
    end

    if p === nothing
        _is_derived_matrix(spec.name) || throw(ArgumentError(
            "Missing required matrix parameter :$(spec.name). Provide a concrete matrix value.",
        ))
        return _derive_matrix(FT, spec.name, ctx, resolved)
    end

    throw(ArgumentError("Internal error: unsupported matrix provider $(typeof(p)) for :$(spec.name)."))
end

# ----------------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------------

"""
    resolve_runtime_parameters(ctx, registry, requirements, ::Type{FT}) -> ModelSpecification

Resolve a `ParamRegistry` to a minimal runtime parameter bundle for the active model.

- Only parameters required by the equations are stored.
- Scalars and vectors are resolved first.
- Matrix parameters are then resolved either from explicit providers or (for a small
  built-in set) via derived defaults.

The return value is a `ModelSpecification` wrapping a `NamedTuple` with keys ordered as:
`(vectors..., matrices..., scalars...)`.
"""
function resolve_runtime_parameters(
    ctx,
    reg::ParamRegistry,
    requirements,
    ::Type{FT},
) where {FT<:AbstractFloat}
    vector_keys = unique(requirements.vectors)
    matrix_keys = unique(requirements.matrices)
    scalar_keys = unique(requirements.scalars)

    # Validate that required parameters exist and have the right declared shape.
    for k in vector_keys
        spec = _require_spec(reg, k)
        spec.shape === :vector || throw(ArgumentError(
            "Parameter :$k is required as a vector but registry declares shape $(spec.shape).",
        ))
    end
    for k in matrix_keys
        spec = _require_spec(reg, k)
        spec.shape === :matrix || throw(ArgumentError(
            "Parameter :$k is required as a matrix but registry declares shape $(spec.shape).",
        ))
    end
    for k in scalar_keys
        spec = _require_spec(reg, k)
        spec.shape === :scalar || throw(ArgumentError(
            "Parameter :$k is required as a scalar but registry declares shape $(spec.shape).",
        ))
    end

    # Determine which scalar/vector keys are needed for resolution.
    needed = Set{Symbol}(vcat(vector_keys, scalar_keys))

    # If a required matrix has no explicit provider and is one of the built-in derived
    # matrices, pull in its declared dependencies.
    for mk in matrix_keys
        spec = _require_spec(reg, mk)
        if spec.provider === nothing
            if haskey(_DERIVED_MATRIX_DEPS, mk)
                for d in _DERIVED_MATRIX_DEPS[mk]
                    push!(needed, d)
                end
            elseif haskey(_DERIVED_MATRIX_DEPS3, mk)
                for d in _DERIVED_MATRIX_DEPS3[mk]
                    push!(needed, d)
                end
            end
        end
    end

    resolved = Dict{Symbol,Any}()

    # Resolve all scalars and vectors in registry order.
    for s in reg.specs
        if s.name in needed
            if s.shape === :scalar
                resolved[s.name] = _resolve_scalar(FT, s, ctx)
            elseif s.shape === :vector
                resolved[s.name] = _resolve_vector(FT, s, ctx)
            end
        end
    end

    # Resolve required matrices (registry order, for stable error messages).
    required_matrices = Set{Symbol}(matrix_keys)
    for s in reg.specs
        if s.shape === :matrix && (s.name in required_matrices)
            resolved[s.name] = _resolve_matrix(FT, s, ctx, resolved)
        end
    end

    # Build minimal runtime NamedTuple in stable order: vectors, matrices, scalars.
    runtime_keys = (vector_keys..., matrix_keys..., scalar_keys...)
    runtime_vals = ntuple(i -> resolved[runtime_keys[i]], length(runtime_keys))

    nt = NamedTuple{Tuple(runtime_keys)}(runtime_vals)
    return ModelSpecification(nt)
end
