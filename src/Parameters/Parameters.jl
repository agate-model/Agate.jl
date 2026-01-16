"""Agate.Parameters

Parameter registries and resolution utilities.

The parameter registry is **CPU-only metadata**. Parameter *resolution* happens during
`construct` and produces a runtime parameter bundle containing only numeric scalars and
arrays. The returned runtime bundle is `Adapt.jl`-compatible, so it can be adapted to
GPU arrays by Oceananigans.

Provider system
---------------
Parameters are declared with an explicit `shape` (`:scalar`, `:vector`, or `:matrix`).
User inputs are normalized into canonical provider wrappers at **boundary time**
(`ParamSpec`, `update_registry`, `extend_registry`, and interactions application).
Runtime resolution never guesses a shape from the value type.

Derived providers may declare dependencies via `deps(provider)::Vector{Symbol}` so the
resolver can compute prerequisite parameters without hardcoded special cases.
"""
module Parameters

using Logging

using Agate.Library.Allometry: AbstractParamDef, resolve_param
using Agate.Library.Allometry:
    allometric_palatability_unimodal_protection,
    palatability_matrix_allometric,
    assimilation_efficiency_matrix_binary

using Agate.Utils.Specifications: PFTSpecification, pft_has, pft_get, ModelSpecification

import Base: show

export ParamSpec, ParamRegistry
export MatrixFn
export scalar_param, vector_param, matrix_param
export parameter_registry, parameter_directory
export resolve_runtime_parameters
export update_registry, extend_registry

# ----------------------------------------------------------------------------
# Shapes
# ----------------------------------------------------------------------------

const _SHAPES = (:scalar, :vector, :matrix)

@inline function _check_shape(shape::Symbol)
    shape in _SHAPES || throw(ArgumentError("Unknown parameter shape $(shape). Expected one of $(_SHAPES)."))
    return nothing
end

# ----------------------------------------------------------------------------
# Provider wrappers (CPU-only)
# ----------------------------------------------------------------------------

"""Abstract supertype for registry providers (CPU-only)."""
abstract type AbstractProvider end

# ---- Scalar providers -------------------------------------------------------

"""Scalar literal provider."""
struct ScalarValue{T} <: AbstractProvider
    value::T
end

"""Scalar provider backed by an allometric definition.

This is only meaningful when used as a *per-PFT* item inside a vector provider.
"""
struct ScalarAllometric{D<:AbstractParamDef} <: AbstractProvider
    def::D
end

"""Scalar provider computed from the construction context: `(ctx) -> value`."""
struct ScalarFromCtx{F} <: AbstractProvider
    f::F
end

# ---- Vector providers -------------------------------------------------------

"""Vector provider storing an explicit full-length vector."""
struct VectorValue{V<:AbstractVector} <: AbstractProvider
    value::V
end

"""Vector provider computed from the construction context: `(ctx) -> vector`."""
struct VectorFromCtx{F} <: AbstractProvider
    f::F
end

"""Vector provider that fills all PFT entries with the same scalar item."""
struct VectorFill <: AbstractProvider
    item::Any
end

"""Vector provider defined by per-group defaults.

The map values are *scalar* providers (applied to all PFTs within that group), and
per-PFT overrides stored in the community specification take precedence.

Accepted inputs for construction:
- `NamedTuple`, e.g. `(Z=1.0, P=0.0)`
- `Dict{Symbol,Any}`
"""
struct VectorGroupMap{M<:AbstractDict{Symbol,Any}} <: AbstractProvider
    map::M
end

# ---- Matrix providers -------------------------------------------------------

"""Matrix provider storing an explicit matrix."""
struct MatrixValue{M<:AbstractMatrix} <: AbstractProvider
    value::M
end


"""\
    MatrixFn(f; deps=Symbol[])

Derived matrix provider with explicit dependencies.

`f` is called during parameter resolution as `f(ctx, deps)` where `deps` is a
`NamedTuple` with fields listed in `deps=...`. Use `deps=[]` for no dependencies.

Use `MatrixFn` to build interaction matrices from other resolved parameters without
hardcoding special cases in the resolver.
"""
struct MatrixFn{F} <: AbstractProvider
    f::F
    deps::Vector{Symbol}
end

MatrixFn(f; deps=Symbol[]) = MatrixFn{typeof(f)}(f, collect(Symbol, deps))

"""Return declared dependencies for a provider."""
deps(::Nothing) = Symbol[]
deps(::AbstractProvider) = Symbol[]
deps(p::MatrixFn) = p.deps

# ----------------------------------------------------------------------------
# Provider normalization (strict-by-default)
# ----------------------------------------------------------------------------

@inline function _normalize_scalar_item(x)
    x === nothing && return nothing

    if x isa ScalarValue || x isa ScalarAllometric
        return x
    elseif x isa AbstractParamDef
        return ScalarAllometric(x)
    elseif x isa Number || x isa Bool
        return ScalarValue(x)
    else
        throw(ArgumentError("Invalid scalar provider item $(typeof(x)). Expected a number/Bool or AbstractParamDef."))
    end
end

function _normalize_vector_provider(x)
    x === nothing && return nothing

    if x isa VectorValue || x isa VectorFromCtx || x isa VectorFill || x isa VectorGroupMap
        return x
    elseif x isa Number || x isa Bool || x isa AbstractParamDef
        return VectorFill(_normalize_scalar_item(x))
    elseif x isa AbstractVector
        return VectorValue(x)
    elseif x isa Function
        return VectorFromCtx(x)
    elseif x isa NamedTuple
        d = Dict{Symbol,Any}()
        for k in propertynames(x)
            d[k] = _normalize_scalar_item(getproperty(x, k))
        end
        return VectorGroupMap(d)
    elseif x isa Dict
        d = Dict{Symbol,Any}()
        for (k, v) in pairs(x)
            k isa Symbol || throw(ArgumentError("Vector group map keys must be Symbol, got $(typeof(k))."))
            d[k] = _normalize_scalar_item(v)
        end
        return VectorGroupMap(d)
    else
        throw(ArgumentError("Invalid vector provider $(typeof(x)). Expected Number/Bool/AbstractParamDef, AbstractVector, NamedTuple/Dict group map, or function."))
    end
end

function _normalize_matrix_provider(x)
    x === nothing && return nothing

    if x isa MatrixValue || x isa MatrixFn
        return x
    elseif x isa AbstractMatrix
        return MatrixValue(x)
    elseif x isa Function
        # Shorthand: treat a bare function as a derived provider with no declared deps.
        # The callable must still support the canonical signature `f(ctx, deps::NamedTuple)`.
        return MatrixFn(x; deps=Symbol[])
    else
        throw(ArgumentError("Invalid matrix provider $(typeof(x)). Expected AbstractMatrix, MatrixFn, or function f(ctx, deps)."))
    end
end

function normalize_provider(shape::Symbol, x)
    _check_shape(shape)

    if shape === :scalar
        x === nothing && return nothing

        if x isa ScalarValue || x isa ScalarFromCtx
            return x
        elseif x isa Function
            return ScalarFromCtx(x)
        elseif x isa Number || x isa Bool
            return ScalarValue(x)
        elseif x isa AbstractParamDef
            throw(ArgumentError("Allometric providers are not valid scalar defaults; use them inside a vector group map or as a vector fill provider."))
        else
            throw(ArgumentError("Invalid scalar provider $(typeof(x)). Expected number/Bool or function."))
        end

    elseif shape === :vector
        return _normalize_vector_provider(x)

    else # :matrix
        return _normalize_matrix_provider(x)
    end
end

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
- `provider`: normalized provider wrapper or `nothing` (meaning required).

Notes
-----
Provider normalization occurs at construction time, so runtime resolution deals only
with canonical provider wrappers.
"""
struct ParamSpec
    name::Symbol
    shape::Symbol
    missing_policy::Symbol
    value_kind::Symbol
    doc::String
    provider::Union{Nothing,AbstractProvider}
end

"""\
    ParamSpec(name, shape, doc, default; missing_policy=:fail, value_kind=:real)

Create a `ParamSpec` and normalize `default` into a canonical provider wrapper.
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

`default` may be a scalar/Bool, a full-length vector, a per-group map (`NamedTuple` or
`Dict{Symbol,Any}`), an allometric definition, or a function `(ctx)->vector`.
"""
vector_param(name::Symbol, doc::AbstractString, default; missing_policy::Symbol=:fail, value_kind::Symbol=:real) =
    ParamSpec(name, :vector, doc, default; missing_policy=missing_policy, value_kind=value_kind)

"""\
    matrix_param(name, doc, default; missing_policy=:fail, value_kind=:real) -> ParamSpec

Create a matrix parameter specification.

`default` may be a concrete matrix, a `MatrixFn(f; deps=[...])`, or a function shorthand `f(ctx, deps)` (normalized to `MatrixFn(f; deps=[])`).
"""
matrix_param(name::Symbol, doc::AbstractString, default; missing_policy::Symbol=:fail, value_kind::Symbol=:real) =
    ParamSpec(name, :matrix, doc, default; missing_policy=missing_policy, value_kind=value_kind)

"""Per-model registry (CPU-only)."""
struct ParamRegistry
    specs::Vector{ParamSpec}
end

"""Return the parameter registry for a factory/model."""
function parameter_registry end

"""Lookup a `ParamSpec` by name."""
function lookup(reg::ParamRegistry, name::Symbol)
    for s in reg.specs
        s.name === name && return s
    end
    return nothing
end

# ----------------------------------------------------------------------------
# Registry updates
# ----------------------------------------------------------------------------

"""\
    update_registry(registry::ParamRegistry; kwargs...) -> ParamRegistry

Return a copy of `registry` with updated parameter providers.

Keys are validated to exist in the registry to catch typos early.
"""
function update_registry(registry::ParamRegistry; kwargs...)
    isempty(kwargs) && return registry

    overrides = (; kwargs...)

    for k in keys(overrides)
        spec = lookup(registry, k)
        spec === nothing && throw(ArgumentError("update_registry: parameter $(k) is not present in this registry"))
    end

    new_specs = Vector{ParamSpec}(undef, length(registry.specs))
    for (i, s) in pairs(registry.specs)
        if hasproperty(overrides, s.name)
            new_default = getproperty(overrides, s.name)
            new_provider = new_default === nothing ? nothing : normalize_provider(s.shape, new_default)
            new_specs[i] = ParamSpec(s.name, s.shape, s.missing_policy, s.value_kind, s.doc, new_provider)
        else
            new_specs[i] = s
        end
    end

    return ParamRegistry(new_specs)
end

"""\
    extend_registry(registry::ParamRegistry, specs::ParamSpec...) -> ParamRegistry

Return a copy of `registry` with additional parameter specifications appended.

To avoid silent mistakes, extending with a key that already exists throws.
"""
function extend_registry(registry::ParamRegistry, specs::ParamSpec...)
    isempty(specs) && return registry

    known = Set(s.name for s in registry.specs)
    for s in specs
        s.name in known && throw(ArgumentError("extend_registry: parameter $(s.name) already exists in registry"))
        push!(known, s.name)
    end

    new_specs = copy(registry.specs)
    append!(new_specs, specs)
    return ParamRegistry(new_specs)
end

# ----------------------------------------------------------------------------
# Directory / pretty printing
# ----------------------------------------------------------------------------

"""Return a lightweight directory of parameters for `factory`."""
function parameter_directory(factory)
    reg = parameter_registry(factory)
    return map(reg.specs) do s
        default_form = isnothing(s.provider) ? :required : typeof(s.provider)
        (name=s.name, shape=s.shape, missing_policy=s.missing_policy, value_kind=s.value_kind, doc=s.doc, default=default_form)
    end
end

@inline _provider_string(p) = isnothing(p) ? "REQUIRED" : string(typeof(p))

function show(io::IO, ::MIME"text/plain", reg::ParamRegistry)
    println(io, "Agate.ParamRegistry with $(length(reg.specs)) parameters")
    for s in reg.specs
        println(io, "  * ", s.name, " (", s.shape, ", ", s.value_kind, ") [", s.missing_policy, "]")
        println(io, "      default: ", _provider_string(s.provider))
        if !isempty(s.doc)
            for line in split(s.doc, '\n')
                println(io, "      ", line)
            end
        end
    end
end

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

function _resolve_scalar_provider(::Type{FT}, spec::ParamSpec, ctx, p) where {FT<:AbstractFloat}
    if p isa ScalarValue
        return coerce_value(FT, spec.value_kind, p.value)    else
        throw(ArgumentError("Internal error: unsupported scalar provider $(typeof(p)) for :$(spec.name)."))
    end
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

    return _resolve_scalar_provider(FT, spec, ctx, p)
end

@inline function _resolve_scalar_item(::Type{FT}, spec::ParamSpec, ctx, idx::Int, p) where {FT<:AbstractFloat}
    p === nothing && return _missing_value(FT, _is_bool(spec))

    if p isa ScalarValue
        return coerce_value(FT, spec.value_kind, p.value)
    elseif p isa ScalarAllometric
        val = resolve_param(FT, p.def, ctx.diameters[idx])
        return coerce_value(FT, spec.value_kind, val)
    else
        throw(ArgumentError("Internal error: unsupported vector item provider $(typeof(p)) for :$(spec.name)."))
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

    out = is_bool ? Vector{Bool}(undef, n) : Vector{FT}(undef, n)

    if p isa VectorValue
        v = p.value
        length(v) == n || throw(ArgumentError("Default for :$(spec.name) must have length $n."))
        @inbounds for i in 1:n
            out[i] = coerce_value(FT, spec.value_kind, v[i])
        end

    elseif p isa VectorFromCtx
        v = p.f(ctx)
        v isa AbstractVector || throw(ArgumentError("Vector provider for :$(spec.name) must return a vector."))
        length(v) == n || throw(ArgumentError("Vector provider for :$(spec.name) must return length $n."))
        @inbounds for i in 1:n
            out[i] = coerce_value(FT, spec.value_kind, v[i])
        end

    elseif p isa VectorFill
        item = p.item
        @inbounds for i in 1:n
            out[i] = _resolve_scalar_item(FT, spec, ctx, i, item)
        end

    elseif p isa VectorGroupMap
        mp = p.map
        @inbounds for i in 1:n
            gsym = ctx.group_symbols[i]
            item = get(mp, gsym, nothing)
            if item === nothing
                _handle_missing!(spec, missing, gsym)
                out[i] = _missing_value(FT, is_bool)
            else
                out[i] = _resolve_scalar_item(FT, spec, ctx, i, item)
            end
        end

    else
        throw(ArgumentError("Internal error: unsupported vector provider $(typeof(p)) for :$(spec.name)."))
    end

    # Apply per-PFT overrides on top of the default.
    @inbounds for i in 1:n
        pft = ctx.pfts[i]
        if pft_has(pft, spec.name)
            raw = pft_get(pft, spec.name)
            if isnothing(raw)
                _handle_missing!(spec, missing, ctx.group_symbols[i])
                out[i] = _missing_value(FT, is_bool)
            else
                item2 = _normalize_scalar_item(raw)
                out[i] = _resolve_scalar_item(FT, spec, ctx, i, item2)
            end
        end
    end

    _emit_missing_warning(spec, missing)
    return out
end

@inline function _deps_namedtuple(deps_syms::Vector{Symbol}, resolved::Dict{Symbol,Any})
    vals = Any[]
    for k in deps_syms
        haskey(resolved, k) || throw(ArgumentError("Matrix provider dependency :$k has not been resolved."))
        push!(vals, resolved[k])
    end
    return NamedTuple{Tuple(deps_syms)}(Tuple(vals))
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

    val = if p isa MatrixValue
        p.value
    elseif p isa MatrixFn
        deps_nt = _deps_namedtuple(p.deps, resolved)
        try
            p.f(ctx, deps_nt)
        catch e
            if e isa MethodError && e.f === p.f
                throw(ArgumentError("MatrixFn for :$(spec.name) must support f(ctx, deps::NamedTuple)."))
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
# Default derived interaction matrices
# ----------------------------------------------------------------------------

@inline function _default_palatability_matrix(ctx, deps)
    return palatability_matrix_allometric(ctx.FT, ctx.diameters;
        can_eat=deps.can_eat,
        optimum_predator_prey_ratio=deps.optimum_predator_prey_ratio,
        specificity=deps.specificity,
        protection=deps.protection,
        palatability_fn=allometric_palatability_unimodal_protection,
    )
end

@inline function _default_assimilation_matrix(ctx, deps)
    FT = eltype(ctx.diameters)
    return assimilation_efficiency_matrix_binary(FT;
        can_eat=deps.can_eat,
        can_be_eaten=deps.can_be_eaten,
        assimilation_efficiency=deps.assimilation_efficiency,
    )
end

"""Return the canonical default palatability matrix provider."""
default_palatability_provider() =
    MatrixFn(_default_palatability_matrix; deps=[:can_eat, :optimum_predator_prey_ratio, :specificity, :protection])

"""Return the canonical default assimilation efficiency matrix provider."""
default_assimilation_provider() =
    MatrixFn(_default_assimilation_matrix; deps=[:can_eat, :can_be_eaten, :assimilation_efficiency])

# ----------------------------------------------------------------------------
# Runtime resolution
# ----------------------------------------------------------------------------

"""Resolve only parameters required by the constructed equations.

Returns a `ModelSpecification` containing only:
- vectors referenced by equations,
- matrices referenced by equations,
- scalars referenced by equations.

Additional parameters needed to build derived providers (via `deps(...)`) are
resolved internally but are not included in the returned runtime bundle.
"""
function resolve_runtime_parameters(
    factory,
    ctx,
    requirements;
    FT::Type{<:AbstractFloat},
    registry=nothing,
)
    reg = (registry === nothing ? parameter_registry(factory) : registry)

    # Build a name->spec index to avoid repeated linear scans.
    specs = Dict{Symbol,ParamSpec}()
    for s in reg.specs
        haskey(specs, s.name) && throw(ArgumentError("Duplicate ParamSpec for :$(s.name) in registry"))
        specs[s.name] = s
    end

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
        spec = get(specs, k, nothing)
        isnothing(spec) && throw(ArgumentError("No ParamSpec found for vector :$k."))
        spec.shape === :vector || throw(ArgumentError("Parameter :$k is required as a vector but registry declares shape $(spec.shape)."))
    end
    for k in matrix_keys
        spec = get(specs, k, nothing)
        isnothing(spec) && throw(ArgumentError("No ParamSpec found for matrix :$k."))
        spec.shape === :matrix || throw(ArgumentError("Parameter :$k is required as a matrix but registry declares shape $(spec.shape)."))
    end
    for k in scalar_keys
        spec = get(specs, k, nothing)
        isnothing(spec) && throw(ArgumentError("No ParamSpec found for scalar :$k."))
        spec.shape === :scalar || throw(ArgumentError("Parameter :$k is required as a scalar but registry declares shape $(spec.shape)."))
    end

    # Compute dependency closure (resolver-only; not part of runtime bundle).
    needed = Set{Symbol}(vcat(vector_keys, matrix_keys, scalar_keys))
    queue = collect(needed)

    while !isempty(queue)
        k = popfirst!(queue)
        spec = get(specs, k, nothing)
        isnothing(spec) && throw(ArgumentError("No ParamSpec found for :$k (required by provider dependencies)."))
        for d in deps(spec.provider)
            if !(d in needed)
                push!(needed, d)
                push!(queue, d)
            end
        end
    end

    resolved = Dict{Symbol,Any}()

    # Resolve all scalars and vectors needed (including dependencies).
    for k in needed
        spec = get(specs, k, nothing)
        spec === nothing && continue
        if spec.shape === :scalar
            resolved[k] = _resolve_scalar(FT, spec, ctx)
        elseif spec.shape === :vector
            resolved[k] = _resolve_vector(FT, spec, ctx)
        end
    end

    # Resolve matrices in dependency order (supports matrix->matrix deps).
    pending = Symbol[]
    for k in needed
        spec = get(specs, k, nothing)
        spec !== nothing && spec.shape === :matrix && push!(pending, k)
    end

    while !isempty(pending)
        progressed = false
        new_pending = Symbol[]

        for k in pending
            spec = specs[k]
            if all(haskey(resolved, d) for d in deps(spec.provider))
                resolved[k] = _resolve_matrix(FT, spec, ctx, resolved)
                progressed = true
            else
                push!(new_pending, k)
            end
        end

        progressed || throw(ArgumentError("Cyclic or unsatisfied matrix provider dependencies: $(new_pending)"))
        pending = new_pending
    end

    # Build minimal runtime NamedTuple in stable order: vectors, matrices, scalars.
    runtime_keys = (vector_keys..., matrix_keys..., scalar_keys...)
    runtime_vals = Any[]

    for k in vector_keys
        push!(runtime_vals, resolved[k])
    end
    for k in matrix_keys
        push!(runtime_vals, resolved[k])
    end
    for k in scalar_keys
        push!(runtime_vals, resolved[k])
    end

    nt = NamedTuple{Tuple(runtime_keys)}(Tuple(runtime_vals))
    return ModelSpecification(nt)
end

end # module Parameters
