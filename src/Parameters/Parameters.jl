"""Agate.Parameters

Parameter registries and resolution utilities for Agate's "new API".

**GPU safety:** the registry is CPU-only metadata. Resolution executes on the CPU and
produces a *runtime* parameter bundle containing only plain numeric scalars and
(Cu)Array-compatible arrays (no Dicts, no functions, no non-isbits structs stored).
"""

module Parameters

using Adapt
using Logging
using Oceananigans

using ..Library.Allometry:
    AbstractParamDef,
    resolve_param,
    allometric_palatability_unimodal_protection,
    palatability_matrix_allometric,
    assimilation_efficiency_matrix_binary

using ..Utils.Specifications: PFTSpecification, pft_has, pft_get, ModelSpecification

export ParamSpec, ParamRegistry
export parameter_registry, parameter_directory
export default_parameter_set, resolve_runtime_parameters
export update_registry, extend_registry
# -----------------------------------------------------------------------------
# Registry types
# -----------------------------------------------------------------------------


"""Single parameter specification.

`default` may be:

* literal scalar/array
* an `AbstractParamDef` (e.g. `AllometricParam`)
* a CPU provider function `(ctx)->value`
* `nothing` meaning "required" (must be overridden when resolved)
"""
struct ParamSpec
    name::Symbol
    # How `nothing` / missing values are treated during resolution.
    #   :fail        -> throw
    #   :zero_warn   -> replace with 0/false and emit a warning
    #   :zero_silent -> replace with 0/false silently
    scope::Symbol
    kind::Symbol
    doc::String
    default::Any
end

"""Convenience constructor.

`scope` controls how `nothing` / missing values are handled during resolution.
"""
ParamSpec(name::Symbol, doc::AbstractString, default; scope::Symbol=:fail, kind::Symbol=:real) =
    ParamSpec(name, scope, kind, String(doc), default)

"""Per-model registry (CPU-only)."""
struct ParamRegistry
    specs::Vector{ParamSpec}
end

"""Return the parameter registry for a factory/model."""
function parameter_registry
end

"""Lookup a `ParamSpec` by name."""
function lookup(reg::ParamRegistry, name::Symbol)
    for s in reg.specs
        s.name === name && return s
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Convenience: non-mutating registry updates
# -----------------------------------------------------------------------------

"""    update_registry(registry::ParamRegistry; kwargs...) -> ParamRegistry

Return a copy of `registry` with updated parameter defaults.

This is intended as the primary user-facing mechanism for overriding parameter
values while keeping a single source of truth (the registry).

Keys are validated to exist in the registry to catch typos early.
"""
function update_registry(registry::ParamRegistry; kwargs...)
    isempty(kwargs) && return registry

    overrides = (; kwargs...)
    for k in keys(overrides)
        lookup(registry, k) === nothing && throw(ArgumentError("update_registry: parameter $(k) is not present in this registry"))
    end

    new_specs = Vector{ParamSpec}(undef, length(registry.specs))
    for (i, s) in pairs(registry.specs)
        if hasproperty(overrides, s.name)
            new_default = getproperty(overrides, s.name)
            new_specs[i] = ParamSpec(s.name, s.scope, s.kind, s.doc, new_default)
        else
            new_specs[i] = s
        end
    end

    return ParamRegistry(new_specs)
end

"""    extend_registry(registry::ParamRegistry, specs::ParamSpec...) -> ParamRegistry

Return a copy of `registry` with additional parameter specifications appended.

This is a lightweight extension hook for model variants and experiments.
"""
function extend_registry(registry::ParamRegistry, specs::ParamSpec...)
    isempty(specs) && return registry
    new_specs = copy(registry.specs)
    append!(new_specs, specs)
    return ParamRegistry(new_specs)
end

# -----------------------------------------------------------------------------
# Value helpers
# -----------------------------------------------------------------------------

@inline function _maybe_ustrip(x)
    # Oceananigans.Units sometimes provides `ustrip` (Unitful-like), but many unit
    # conveniences are simple numeric scale factors (day, hour, second).
    if isdefined(Oceananigans.Units, :ustrip)
        try
            return Oceananigans.Units.ustrip(x)
        catch
            return x
        end
    end
    return x
end

@inline strip_units(x::Number) = _maybe_ustrip(x)
@inline strip_units(x::AbstractArray) = strip_units.(x)
@inline strip_units(x) = x

# -----------------------------------------------------------------------------
# Directory / introspection
# -----------------------------------------------------------------------------

"""Return a lightweight directory of parameters for `factory`.

This is CPU-only introspection: suitable for printing/saving/diffing.
"""
function parameter_directory(factory)
    reg = parameter_registry(factory)
    return map(reg.specs) do s
        provider = s.default
        default_form = isnothing(provider) ? :required : typeof(provider)
        (name=s.name, scope=s.scope, kind=s.kind, doc=s.doc, default=default_form)
    end
end

# -----------------------------------------------------------------------------
# Resolution
# -----------------------------------------------------------------------------


@inline _is_bool(spec::ParamSpec) = spec.kind === :bool

@inline _missing_value(::Type{FT}, is_bool::Bool) where {FT<:AbstractFloat} = is_bool ? false : zero(FT)

@inline function _handle_missing!(spec::ParamSpec, missing::Vector{Symbol}, gsym::Symbol)
    if spec.scope === :fail
        throw(ArgumentError("Missing value for parameter :$(spec.name) (group $(gsym))."))
    elseif spec.scope === :zero_warn
        push!(missing, gsym)
    end
    return nothing
end

@inline function _emit_missing_warning(spec::ParamSpec, missing::Vector{Symbol})
    if spec.scope === :zero_warn && !isempty(missing)
        # Warn once per resolution call; include unique group symbols.
        @warn "Parameter :$(spec.name) missing for groups $(unique(missing)); replacing with 0/false." maxlog=1
    end
    return nothing
end

"""Coerce a scalar-like value to `FT` after unit stripping."""
@inline function _to_FT(::Type{FT}, x) where {FT<:AbstractFloat}
    y = strip_units(x)
    return y isa Bool ? y : FT(y)
end

"""Resolve a single scalar parameter (CPU)."""
function _resolve_scalar(::Type{FT}, spec::ParamSpec, ctx) where {FT<:AbstractFloat}
    provider = spec.default

    if isnothing(provider)
        if spec.scope === :fail
            throw(ArgumentError("Missing required scalar parameter :$(spec.name)."))
        elseif spec.scope === :zero_warn
            @warn "Missing scalar parameter :$(spec.name); replacing with 0." maxlog=1
        end
        return _missing_value(FT, _is_bool(spec))
    end

    val = provider isa Function ? provider(ctx) : provider
    return _to_FT(FT, val)
end

"""Resolve a single element for a vector parameter."""
@inline function _resolve_vector_element(::Type{FT}, provider, diameter) where {FT<:AbstractFloat}
    if provider isa AbstractParamDef
        # resolve_param already handles casting to FT for numeric outputs.
        return _to_FT(FT, resolve_param(FT, provider, diameter))
    elseif provider isa Bool
        return provider
    else
        return _to_FT(FT, provider)
    end
end

@inline function _group_provider(provider, gsym::Symbol)
    if provider isa NamedTuple
        return hasproperty(provider, gsym) ? getproperty(provider, gsym) : nothing
    elseif provider isa Dict
        return get(provider, gsym, nothing)
    else
        return provider
    end
end

"""Resolve a vector parameter over the full plankton community.

Precedence:
1. explicit per-PFT overrides stored in `ctx.pfts[i]`
2. registry default (either a full vector or a per-group mapping)
3. fallback: zero/false (for groups missing from a mapping)
"""
function _resolve_vector(::Type{FT}, spec::ParamSpec, ctx) where {FT<:AbstractFloat}
    n = ctx.n_total
    is_bool = _is_bool(spec)
    missing = Symbol[]

    default_provider = spec.default
    if isnothing(default_provider)
        if spec.scope === :fail
            throw(ArgumentError("Missing required vector parameter :$(spec.name)."))
        elseif spec.scope === :zero_warn
            @warn "Missing vector parameter :$(spec.name); replacing with 0/false." maxlog=1
        end
        return is_bool ? fill(false, n) : fill(zero(FT), n)
    end

    # Full-vector defaults are allowed in the registry (useful for explicit
    # overrides via `update_registry`). Apply per-PFT overrides on top.
    if default_provider isa AbstractVector
        length(default_provider) == n || throw(ArgumentError("Default for :$(spec.name) must have length $n."))
        out = is_bool ? Bool.(default_provider) : _to_FT.(Ref(FT), default_provider)
        @inbounds for i in 1:n
            pft = ctx.pfts[i]
            if pft_has(pft, spec.name)
                v = pft_get(pft, spec.name)
                if isnothing(v)
                    _handle_missing!(spec, missing, ctx.group_symbols[i])
                    out[i] = _missing_value(FT, is_bool)
                else
                    v = v isa Function ? v(ctx) : v
                    out[i] = _resolve_vector_element(FT, v, ctx.diameters[i])
                end
            end
        end
        _emit_missing_warning(spec, missing)
        return out
    end

    out = is_bool ? Vector{Bool}(undef, n) : Vector{FT}(undef, n)
    @inbounds for i in 1:n
        pft = ctx.pfts[i]
        if pft_has(pft, spec.name)
            v = pft_get(pft, spec.name)
            if isnothing(v)
                _handle_missing!(spec, missing, ctx.group_symbols[i])
                out[i] = _missing_value(FT, is_bool)
            else
                v = v isa Function ? v(ctx) : v
                out[i] = _resolve_vector_element(FT, v, ctx.diameters[i])
            end
            continue
        end

        gsym = ctx.group_symbols[i]
        p = _group_provider(default_provider, gsym)
        if isnothing(p)
            _handle_missing!(spec, missing, gsym)
            out[i] = _missing_value(FT, is_bool)
        else
            p = p isa Function ? p(ctx) : p
            out[i] = _resolve_vector_element(FT, p, ctx.diameters[i])
        end
    end

    _emit_missing_warning(spec, missing)

    return out
end

"""Default palatability matrix provider (CPU).

`palatability_fn` is the *pairwise* palatability kernel used to build the full matrix.
"""
function default_palatability_matrix(ctx, vectors; palatability_fn=allometric_palatability_unimodal_protection)
    return palatability_matrix_allometric(ctx.FT, ctx.diameters;
        can_eat=vectors[:can_eat],
        optimum_predator_prey_ratio=vectors[:optimum_predator_prey_ratio],
        specificity=vectors[:specificity],
        protection=vectors[:protection],
        palatability_fn=palatability_fn,
    )
end

"""Default assimilation efficiency matrix provider (CPU)."""
function default_assimilation_matrix(ctx, vectors)
    FT = eltype(ctx.diameters)
    can_eat   = vectors[:can_eat]
    can_be    = vectors[:can_be_eaten]
    assim_eff = vectors[:assimilation_efficiency]
    return assimilation_efficiency_matrix_binary(FT; can_eat, can_be_eaten=can_be, assimilation_efficiency=assim_eff)
end

function _resolve_matrix(::Type{FT}, spec::ParamSpec, ctx, vectors;
                         palatability_fn=allometric_palatability_unimodal_protection) where {FT<:AbstractFloat}
    n = ctx.n_total
    provider = spec.default
    if isnothing(provider)
        if spec.scope === :fail
            throw(ArgumentError("Missing required matrix parameter :$(spec.name)."))
        elseif spec.scope === :zero_warn
            @warn "Missing matrix parameter :$(spec.name); replacing with zeros." maxlog=1
        end
        return zeros(FT, n, n)
    end

    # Compute providers execute on CPU during resolution only.
    if provider === default_palatability_matrix
        val = default_palatability_matrix(ctx, vectors; palatability_fn)
    elseif provider === default_assimilation_matrix
        val = default_assimilation_matrix(ctx, vectors)
    elseif provider isa Function
        val = provider(ctx)
    else
        val = provider
    end

    (val isa AbstractMatrix) || throw(ArgumentError("Matrix parameter :$(spec.name) must be a matrix."))
    size(val, 1) == n && size(val, 2) == n || throw(ArgumentError("Matrix :$(spec.name) must be size ($n,$n)."))

    return FT.(strip_units.(val))
end

"""Generate a complete resolved default set for `factory` (CPU arrays/scalars).

This is meant for printing/saving/diffing; it resolves *all* registry entries.
"""
function default_parameter_set(factory, ctx;
    FT::Type{<:AbstractFloat},
    palatability_fn=allometric_palatability_unimodal_protection,
)
    reg = parameter_registry(factory)

    # Classify each spec for printing/saving by looking at its selected provider.
    # Runtime construction does *not* use this; it uses `requirements`.
    function infer_rank(spec::ParamSpec)
        provider = spec.default
        if spec.name === :palatability_matrix || spec.name === :assimilation_matrix
            return :matrix
        elseif provider === default_palatability_matrix || provider === default_assimilation_matrix
            return :matrix
        elseif provider isa AbstractMatrix
            return :matrix
        elseif occursin("matrix", String(spec.name)) && provider isa Function
            return :matrix
        elseif provider isa AbstractVector || provider isa NamedTuple || provider isa Dict || provider isa AbstractParamDef
            return :vector
        else
            return :scalar
        end
    end

    scalars = Dict{Symbol,Any}()
    vectors = Dict{Symbol,Any}()
    matrices = Dict{Symbol,Any}()

    ranks = Dict{Symbol,Symbol}()
    for spec in reg.specs
        r = infer_rank(spec)
        ranks[spec.name] = r
        if r === :scalar
            scalars[spec.name] = _resolve_scalar(FT, spec, ctx)
        elseif r === :vector
            vectors[spec.name] = _resolve_vector(FT, spec, ctx)
        end
    end

    for spec in reg.specs
        if ranks[spec.name] === :matrix
            matrices[spec.name] = _resolve_matrix(FT, spec, ctx, vectors; palatability_fn)
        end
    end

    # Return a NamedTuple in registry order (stable for printing/saving/diffing).
    out_keys = Symbol[]
    out_vals = Any[]
    for spec in reg.specs
        r = ranks[spec.name]
        push!(out_keys, spec.name)
        if r === :scalar
            push!(out_vals, scalars[spec.name])
        elseif r === :vector
            push!(out_vals, vectors[spec.name])
        else
            push!(out_vals, matrices[spec.name])
        end
    end

    return NamedTuple{Tuple(out_keys)}(Tuple(out_vals))
end

"""Resolve only parameters required by the constructed equations.

Returns a `ModelSpecification` with a minimal set of runtime parameters.
"""
function resolve_runtime_parameters(
    factory,
    ctx,
    requirements;
    FT::Type{<:AbstractFloat},
    palatability_fn=allometric_palatability_unimodal_protection,
    registry=nothing,
)
    reg = (registry === nothing ? parameter_registry(factory) : registry)

    # What the equations actually need.
    vector_keys = unique(requirements.vectors)
    matrix_keys = requirements.matrices
    scalar_keys = requirements.scalars

    # Resolve only vectors needed by equations + any dependencies for core defaults.
    dep_vecs = Symbol[]
    if (:palatability_matrix in matrix_keys)
        spec = lookup(reg, :palatability_matrix)
        if !isnothing(spec) && spec.default === default_palatability_matrix
            append!(dep_vecs, (:can_eat, :optimum_predator_prey_ratio, :specificity, :protection))
        end
    end
    if (:assimilation_matrix in matrix_keys)
        spec = lookup(reg, :assimilation_matrix)
        if !isnothing(spec) && spec.default === default_assimilation_matrix
            append!(dep_vecs, (:can_eat, :can_be_eaten, :assimilation_efficiency))
        end
    end

    needed_vectors = unique(vcat(vector_keys, dep_vecs))

    # Resolve scalars.
    scalar_vals = Dict{Symbol,Any}()
    for k in scalar_keys
        spec = lookup(reg, k)
        isnothing(spec) && throw(ArgumentError("No ParamSpec found for scalar :$k."))
        scalar_vals[k] = _resolve_scalar(FT, spec, ctx)
    end

    # Resolve vectors (including dependencies).
    vector_vals = Dict{Symbol,Any}()
    for k in needed_vectors
        spec = lookup(reg, k)
        isnothing(spec) && throw(ArgumentError("No ParamSpec found for vector :$k."))
        vector_vals[k] = _resolve_vector(FT, spec, ctx)
    end

    # Resolve matrices.
    matrix_vals = Dict{Symbol,Any}()
    for k in matrix_keys
        spec = lookup(reg, k)
        isnothing(spec) && throw(ArgumentError("No ParamSpec found for matrix :$k."))
        matrix_vals[k] = _resolve_matrix(FT, spec, ctx, vector_vals; palatability_fn)
    end

    # Build minimal runtime NamedTuple in stable order: vectors, matrices, scalars.
    # (Order only matters for printing; tracer binding uses propertynames.)
    runtime_keys = (vector_keys..., matrix_keys..., scalar_keys...)
    runtime_vals = Any[]
    for k in vector_keys
        push!(runtime_vals, vector_vals[k])
    end
    for k in matrix_keys
        push!(runtime_vals, matrix_vals[k])
    end
    for k in scalar_keys
        push!(runtime_vals, scalar_vals[k])
    end

    nt = NamedTuple{Tuple(runtime_keys)}(Tuple(runtime_vals))
    params = ModelSpecification(nt)

    return params
end


end # module Parameters
