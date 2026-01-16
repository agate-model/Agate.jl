# This file is included into the `Agate.Parameters` module.

# ----------------------------------------------------------------------------
# Resolution helpers
# ----------------------------------------------------------------------------

@inline _is_bool(spec::ParamSpec) = spec.value_kind === :bool

"""Convert a user value to the runtime storage type for a given `value_kind`."""
@inline function coerce_value(::Type{FT}, value_kind::Symbol, x) where {FT<:AbstractFloat}
    value_kind === :bool && return Bool(x)
    return FT(x)
end

@inline function _missing_value(::Type{FT}, is_bool::Bool) where {FT<:AbstractFloat}
    return is_bool ? false : zero(FT)
end

@inline function _handle_missing!(spec::ParamSpec, missing_groups::Vector{Symbol}, group::Symbol)
    if spec.missing_policy === :fail
        throw(ArgumentError("Missing :$(spec.name) for group :$(group). Provide an override or change missing_policy."))
    elseif spec.missing_policy === :zero_warn
        push!(missing_groups, group)
    end
    return nothing
end

@inline function _emit_missing_warning(spec::ParamSpec, missing_groups::Vector{Symbol})
    if spec.missing_policy === :zero_warn && !isempty(missing_groups)
        uniq = unique(missing_groups)
        @warn "Missing entries for vector :$(spec.name); replacing with 0/false for groups $(uniq)." maxlog=1
    end
    return nothing
end

function _resolve_scalar(::Type{FT}, spec::ParamSpec, ctx) where {FT<:AbstractFloat}
    spec.shape === :scalar || throw(ArgumentError("Internal error: :$(spec.name) is not a scalar spec."))

    p = spec.provider
    if isnothing(p)
        if spec.missing_policy === :fail
            throw(ArgumentError("Missing required scalar parameter :$(spec.name)."))
        elseif spec.missing_policy === :zero_warn
            @warn "Missing scalar parameter :$(spec.name); replacing with 0/false." maxlog=1
        end
        return _missing_value(FT, _is_bool(spec))
    end

    if p isa Number || p isa Bool
        return coerce_value(FT, spec.value_kind, p)
    end

    throw(ArgumentError("Internal error: unsupported scalar provider $(typeof(p)) for :$(spec.name)."))
end

@inline function _resolve_scalar_item(::Type{FT}, spec::ParamSpec, ctx, idx::Int, p) where {FT<:AbstractFloat}
    p === nothing && return _missing_value(FT, _is_bool(spec))

    if p isa Number || p isa Bool
        return coerce_value(FT, spec.value_kind, p)
    elseif p isa AbstractParamDef
        val = resolve_param(FT, p, ctx.diameters[idx])
        return coerce_value(FT, spec.value_kind, val)
    end

    throw(ArgumentError("Internal error: unsupported vector item provider $(typeof(p)) for :$(spec.name)."))
end

@inline function _default_vector_item(spec::ParamSpec, provider, ctx, i::Int)
    if provider isa AbstractVector
        return true, provider[i]
    elseif provider isa Number || provider isa Bool || provider isa AbstractParamDef
        return true, provider
    elseif provider isa VectorGroupMap
        gsym = ctx.group_symbols[i]
        ks = provider.keys
        its = provider.items
        for j in eachindex(ks)
            if ks[j] === gsym
                return true, its[j]
            end
        end
        return false, nothing
    else
        throw(ArgumentError("Internal error: unsupported vector provider $(typeof(provider)) for :$(spec.name)."))
    end
end

"""Resolve a vector parameter over the full plankton community.

Precedence:
1. explicit per-PFT overrides stored in `ctx.pfts[i]`
2. registry default provider
3. fallback: 0/false for missing groups (according to `missing_policy`)
"""
function _resolve_vector(::Type{FT}, spec::ParamSpec, ctx) where {FT<:AbstractFloat}
    spec.shape === :vector || throw(ArgumentError("Internal error: :$(spec.name) is not a vector spec."))

    n = ctx.n_total
    is_bool = _is_bool(spec)
    missing = Symbol[]

    p = spec.provider
    if isnothing(p)
        if spec.missing_policy === :fail
            throw(ArgumentError("Missing required vector parameter :$(spec.name)."))
        elseif spec.missing_policy === :zero_warn
            @warn "Missing vector parameter :$(spec.name); replacing with 0/false." maxlog=1
        end
        return is_bool ? fill(false, n) : fill(zero(FT), n)
    end

    if p isa AbstractVector
        length(p) == n || throw(ArgumentError("Default for :$(spec.name) must have length $n."))
    end

    out = is_bool ? Vector{Bool}(undef, n) : Vector{FT}(undef, n)

    @inbounds for i in 1:n
        pft = ctx.pfts[i]

        raw = if pft_has(pft, spec.name)
            pft_get(pft, spec.name)
        else
            found, item = _default_vector_item(spec, p, ctx, i)
            if !found
                _handle_missing!(spec, missing, ctx.group_symbols[i])
                nothing
            else
                item
            end
        end

        if isnothing(raw)
            _handle_missing!(spec, missing, ctx.group_symbols[i])
            out[i] = _missing_value(FT, is_bool)
        else
            item2 = pft_has(pft, spec.name) ? _normalize_scalar_item(raw) : raw
            out[i] = _resolve_scalar_item(FT, spec, ctx, i, item2)
        end
    end

    _emit_missing_warning(spec, missing)
    return out
end

function _resolve_matrix(::Type{FT}, spec::ParamSpec, ctx, resolved::Dict{Symbol,Any}) where {FT<:AbstractFloat}
    spec.shape === :matrix || throw(ArgumentError("Internal error: :$(spec.name) is not a matrix spec."))

    n = ctx.n_total
    p = spec.provider

    if isnothing(p)
        if spec.missing_policy === :fail
            throw(ArgumentError("Missing required matrix parameter :$(spec.name)."))
        elseif spec.missing_policy === :zero_warn
            @warn "Missing matrix parameter :$(spec.name); replacing with zeros." maxlog=1
        end
        return zeros(FT, n, n)
    end

    val = if p isa AbstractMatrix
        p
    elseif p isa MatrixFn
        ds = deps(p)
        depvals = ntuple(i -> begin
            k = ds[i]
            haskey(resolved, k) || throw(ArgumentError("Matrix provider dependency :$k has not been resolved."))
            resolved[k]
        end, length(ds))
        try
            p.f(ctx, depvals)
        catch e
            if e isa MethodError && e.f === p.f
                throw(ArgumentError("MatrixFn for :$(spec.name) must support f(ctx, depvals::Tuple)."))
            end
            rethrow()
        end
    else
        throw(ArgumentError("Internal error: unsupported matrix provider $(typeof(p)) for :$(spec.name)."))
    end

    (val isa AbstractMatrix) || throw(ArgumentError("Matrix parameter :$(spec.name) must be a matrix."))
    size(val, 1) == n && size(val, 2) == n || throw(ArgumentError("Matrix :$(spec.name) must be size ($n,$n)."))

    return spec.value_kind === :bool ? Bool.(val) : FT.(val)
end

# ----------------------------------------------------------------------------
# Runtime resolution
# ----------------------------------------------------------------------------

"""Resolve only parameters required by the constructed equations.

Returns a `ModelSpecification` containing only:
- vectors referenced by equations,
- matrices referenced by equations,
- scalars referenced by equations.

Additional parameters needed to build derived providers (via `deps(...)`) are resolved
internally but are not included in the returned runtime bundle.
"""
function resolve_runtime_parameters(
    factory,
    ctx,
    requirements;
    FT::Type{<:AbstractFloat},
    registry=nothing,
)
    reg = (registry === nothing ? parameter_registry(factory) : registry)

    # Guardrail: a parameter key must not be required in multiple shapes.
    all_keys = vcat(requirements.scalars, requirements.vectors, requirements.matrices)
    if length(all_keys) != length(unique(all_keys))
        dup = [k for k in unique(all_keys) if count(==(k), all_keys) > 1]
        throw(ArgumentError("Parameter keys appear in multiple shapes (scalar/vector/matrix): $(dup)"))
    end

    vector_keys = unique(requirements.vectors)
    matrix_keys = unique(requirements.matrices)
    scalar_keys = unique(requirements.scalars)

    # Validate declared shapes.
    for k in vector_keys
        spec = _require_spec(reg, k)
        spec.shape === :vector || throw(ArgumentError("Parameter :$k is required as a vector but registry declares shape $(spec.shape)."))
    end
    for k in matrix_keys
        spec = _require_spec(reg, k)
        spec.shape === :matrix || throw(ArgumentError("Parameter :$k is required as a matrix but registry declares shape $(spec.shape)."))
    end
    for k in scalar_keys
        spec = _require_spec(reg, k)
        spec.shape === :scalar || throw(ArgumentError("Parameter :$k is required as a scalar but registry declares shape $(spec.shape)."))
    end

    # Compute dependency closure (resolver-only; not part of runtime bundle).
    needed = Set{Symbol}(vcat(vector_keys, matrix_keys, scalar_keys))
    stack = collect(needed)

    while !isempty(stack)
        k = pop!(stack)
        spec = _require_spec(reg, k)
        for d in deps(spec.provider)
            dep_spec = _require_spec(reg, d)
            if spec.shape === :matrix && dep_spec.shape === :matrix
                throw(ArgumentError(
                    "Matrix parameter :$k declares matrix dependency :$d. MatrixFn providers may only depend on scalar/vector parameters.",
                ))
            end
            if !(d in needed)
                push!(needed, d)
                push!(stack, d)
            end
        end
    end

    resolved = Dict{Symbol,Any}()

    # Resolve all scalars and vectors needed (including dependencies) in registry order.
    for s in reg.specs
        if s.name in needed
            if s.shape === :scalar
                resolved[s.name] = _resolve_scalar(FT, s, ctx)
            elseif s.shape === :vector
                resolved[s.name] = _resolve_vector(FT, s, ctx)
            end
        end
    end

    # Resolve matrices. Matrix providers are restricted to scalar/vector dependencies,
    # so a single pass in registry order is sufficient.
    for s in reg.specs
        if s.shape === :matrix && (s.name in needed)
            resolved[s.name] = _resolve_matrix(FT, s, ctx, resolved)
        end
    end

    # Build minimal runtime NamedTuple in stable order: vectors, matrices, scalars.
    runtime_keys = (vector_keys..., matrix_keys..., scalar_keys...)
    runtime_vals = ntuple(i -> resolved[runtime_keys[i]], length(runtime_keys))

    nt = NamedTuple{Tuple(runtime_keys)}(runtime_vals)
    return ModelSpecification(nt)
end
