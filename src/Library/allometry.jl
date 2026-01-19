module Allometry

"""Utilities for size-dependent traits and interaction matrices."""

export AbstractParamDef, ConstantParam, AllometricParam
export PowerLaw, resolve_param, cast_paramdef

export PalatabilityPreyParameters, PalatabilityPredatorParameters
export allometric_scaling_power
export allometric_palatability_unimodal, allometric_palatability_unimodal_protection
export palatability_matrix_allometric, assimilation_efficiency_matrix_binary
export palatability_matrix_allometric_axes, assimilation_efficiency_matrix_binary_axes
export build_palatability_matrix, build_assimilation_matrix

# -----------------------------------------------------------------------------
# Explicit parameter definitions
# -----------------------------------------------------------------------------

"""Abstract supertype for explicit parameter definitions stored in `PFTSpecification`."""
abstract type AbstractParamDef end

"""A parameter that is constant across sizes."""
struct ConstantParam{T} <: AbstractParamDef
    value::T
end

ConstantParam(x) = ConstantParam{typeof(x)}(x)

"""A parameter that is computed from an explicit model/callable and coefficient bundle."""
struct AllometricParam{F,C} <: AbstractParamDef
    model::F
    coeffs::C
end

"""Construct an `AllometricParam` from a model and keyword coefficients."""
AllometricParam(model; kwargs...) = AllometricParam(model, (; kwargs...))

"""Common allometry model: power-law scaling with spherical volume."""
struct PowerLaw end

"""Evaluate a `PowerLaw` model.

Expected coefficient names:
- `prefactor` and `exponent` (recommended), or
- `scale` and `exponent` (alias), or
- `a` and `b` (legacy).
"""
@inline function (m::PowerLaw)(coeffs::NamedTuple, diameter)
    if hasproperty(coeffs, :prefactor)
        a = getproperty(coeffs, :prefactor)
    elseif hasproperty(coeffs, :scale)
        a = getproperty(coeffs, :scale)
    elseif hasproperty(coeffs, :a)
        a = getproperty(coeffs, :a)
    else
        throw(ArgumentError("PowerLaw requires coefficient `prefactor` (or `scale`/`a`)"))
    end

    if hasproperty(coeffs, :exponent)
        b = getproperty(coeffs, :exponent)
    elseif hasproperty(coeffs, :b)
        b = getproperty(coeffs, :b)
    else
        throw(ArgumentError("PowerLaw requires coefficient `exponent` (or `b`)"))
    end

    # By construction we keep coefficients and diameter the same floating type (FT),
    # so this call never mixes Float32/Float64 (important for GPU use).
    return allometric_scaling_power(a, b, diameter)
end

"""Resolve a parameter definition at a given diameter.

This is the one function the constructor uses when building runtime parameter vectors.
"""
@inline resolve_param(::Type{FT}, x, diameter) where {FT<:AbstractFloat} = FT(x)

@inline resolve_param(::Type{FT}, x::Bool, diameter) where {FT<:AbstractFloat} = x

@inline resolve_param(::Type{FT}, p::ConstantParam, diameter) where {FT<:AbstractFloat} = FT(p.value)

@inline function resolve_param(::Type{FT}, p::AllometricParam, diameter) where {FT<:AbstractFloat}
    # Coefficients often come from literal numbers (Float64). Convert them to FT so we
    # never mix Float32/Float64 in the underlying allometric calls.
    coeffs = map(v -> v isa Number ? FT(v) : v, p.coeffs)
    return FT(p.model(coeffs, FT(diameter)))
end

"""Cast numeric entries inside a parameter definition to `FT`."""
@inline function cast_paramdef(::Type{FT}, p::ConstantParam) where {FT<:AbstractFloat}
    return ConstantParam(FT(p.value))
end

@inline function cast_paramdef(::Type{FT}, p::AllometricParam) where {FT<:AbstractFloat}
    coeffs = p.coeffs
    if coeffs isa NamedTuple
        coeffs = map(v -> v isa Number ? FT(v) : v, coeffs)
    end
    return AllometricParam(p.model, coeffs)
end

# -----------------------------------------------------------------------------
# Palatability + assimilation matrix utilities
# -----------------------------------------------------------------------------

struct PalatabilityPreyParameters{FT<:AbstractFloat}
    diameter::FT
    protection::FT
end

struct PalatabilityPredatorParameters{FT<:AbstractFloat}
    can_eat::Bool
    diameter::FT
    optimum_predator_prey_ratio::FT
    specificity::FT
end

"""Power-law scaling function using spherical volume."""
@inline function allometric_scaling_power(a::FT, b::FT, diameter::FT) where {FT<:AbstractFloat}
    r = diameter / FT(2)
    volume = (FT(4) / FT(3)) * FT(π) * r^FT(3)
    return a * volume^b
end

"""Unimodal palatability (no protection)."""
@inline function allometric_palatability_unimodal(
    prey::PalatabilityPreyParameters{FT}, predator::PalatabilityPredatorParameters{FT}
) where {FT<:AbstractFloat}
    predator.can_eat || return zero(FT)
    ratio = predator.diameter / prey.diameter
    width = one(FT) + (ratio - predator.optimum_predator_prey_ratio)^FT(2)
    return one(FT) / width^predator.specificity
end

"""Unimodal palatability with prey protection.

Protection η reduces palatability as `(1 - η)` (so η=0 means no protection).
"""
 function allometric_palatability_unimodal_protection(
    prey::PalatabilityPreyParameters{FT}, predator::PalatabilityPredatorParameters{FT}
) where {FT<:AbstractFloat}
    base = allometric_palatability_unimodal(prey, predator)
    return base * (one(FT) - prey.protection)
end

"""Build an allometric palatability matrix `M[pred, prey]`."""
function palatability_matrix_allometric(
    ::Type{FT},
    diameters::AbstractVector{FT};
    can_eat::AbstractVector{Bool},
    can_be_eaten::AbstractVector{Bool},
    optimum_predator_prey_ratio::AbstractVector{FT},
    specificity::AbstractVector{FT},
    protection::AbstractVector{FT},
    palatability_fn=allometric_palatability_unimodal_protection,
) where {FT<:AbstractFloat}
    n = length(diameters)
    M = zeros(FT, n, n)

    @inbounds for pred in 1:n
        if !can_eat[pred]
            continue
        end

        predator = PalatabilityPredatorParameters{FT}(
            true,
            diameters[pred],
            optimum_predator_prey_ratio[pred],
            specificity[pred],
        )

        for prey in 1:n
            if !can_be_eaten[prey]
                M[pred, prey] = zero(FT)
                continue
            end
            prey_params = PalatabilityPreyParameters{FT}(diameters[prey], protection[prey])
            M[pred, prey] = palatability_fn(prey_params, predator)
        end
    end

    return M
end

"""Build a role-aware allometric palatability matrix `M[consumer, prey]`.

This computes only the consumer-by-prey block specified by `consumer_indices`
and `prey_indices`, avoiding allocation of an intermediate full square matrix.
"""
function palatability_matrix_allometric_axes(
    ::Type{FT},
    diameters::AbstractVector{FT};
    can_eat::AbstractVector{Bool},
    can_be_eaten::AbstractVector{Bool},
    optimum_predator_prey_ratio::AbstractVector{FT},
    specificity::AbstractVector{FT},
    protection::AbstractVector{FT},
    consumer_indices::AbstractVector{<:Integer},
    prey_indices::AbstractVector{<:Integer},
    palatability_fn=allometric_palatability_unimodal_protection,
) where {FT<:AbstractFloat}
    nr = length(consumer_indices)
    nc = length(prey_indices)
    M = zeros(FT, nr, nc)

    @inbounds for (ii, pred) in pairs(consumer_indices)
        if !can_eat[pred]
            continue
        end

        predator = PalatabilityPredatorParameters{FT}(
            true,
            diameters[pred],
            optimum_predator_prey_ratio[pred],
            specificity[pred],
        )

        for (jj, prey) in pairs(prey_indices)
            if !can_be_eaten[prey]
                M[ii, jj] = zero(FT)
                continue
            end
            prey_params = PalatabilityPreyParameters{FT}(diameters[prey], protection[prey])
            M[ii, jj] = palatability_fn(prey_params, predator)
        end
    end

    return M
end

"""Build a binary assimilation-efficiency matrix `β[pred, prey]`."""
function assimilation_efficiency_matrix_binary(
    ::Type{FT};
    can_eat::AbstractVector{Bool},
    can_be_eaten::AbstractVector{Bool},
    assimilation_efficiency::AbstractVector{FT},
) where {FT<:AbstractFloat}
    n = length(can_eat)
    M = zeros(FT, n, n)
    @inbounds for pred in 1:n
        if !can_eat[pred]
            continue
        end
        β = assimilation_efficiency[pred]
        for prey in 1:n
            M[pred, prey] = can_be_eaten[prey] ? β : zero(FT)
        end
    end
    return M
end

"""Build a role-aware binary assimilation-efficiency matrix `β[consumer, prey]`.

Only the consumer-by-prey block specified by `consumer_indices` and
`prey_indices` is computed.
"""
function assimilation_efficiency_matrix_binary_axes(
    ::Type{FT};
    can_eat::AbstractVector{Bool},
    can_be_eaten::AbstractVector{Bool},
    assimilation_efficiency::AbstractVector{FT},
    consumer_indices::AbstractVector{<:Integer},
    prey_indices::AbstractVector{<:Integer},
) where {FT<:AbstractFloat}
    nr = length(consumer_indices)
    nc = length(prey_indices)
    M = zeros(FT, nr, nc)

    @inbounds for (ii, pred) in pairs(consumer_indices)
        if !can_eat[pred]
            continue
        end
        β = assimilation_efficiency[pred]
        for (jj, prey) in pairs(prey_indices)
            M[ii, jj] = can_be_eaten[prey] ? β : zero(FT)
        end
    end

    return M
end

"""Convenience overload that infers `FT` from `assimilation_efficiency`.

This avoids call sites needing to thread `FT` explicitly while still preserving
GPU compatibility (no mixed float types).
"""
function assimilation_efficiency_matrix_binary(
    can_eat::AbstractVector{Bool},
    can_be_eaten::AbstractVector{Bool},
    assimilation_efficiency::AbstractVector{FT},
) where {FT<:AbstractFloat}
    return assimilation_efficiency_matrix_binary(
        FT;
        can_eat=can_eat,
        can_be_eaten=can_be_eaten,
        assimilation_efficiency=assimilation_efficiency,
    )
end

# -----------------------------------------------------------------------------
# High-level interaction matrix builders
# -----------------------------------------------------------------------------

@inline function _trait_get(pft_data, key::Symbol, default)
    if !hasproperty(pft_data, key)
        return default
    end
    v = getproperty(pft_data, key)
    return v === nothing ? default : v
end

@inline function _resolve_trait(::Type{FT}, v, diameter::FT) where {FT<:AbstractFloat}
    return resolve_param(FT, v, diameter)
end

"""
    build_palatability_matrix(FT, pft_data, diameters;
                              palatability_fn=allometric_palatability_unimodal_protection)

Build the default palatability matrix `M[pred, prey]` from PFT trait definitions.

`pft_data` must be an indexable collection (length `n`) where each element supports
`hasproperty`/`getproperty` for the trait keys used by the default rule:

- `can_eat::Bool` (default `false`)
- `can_be_eaten::Bool` (default `true`)
- `optimum_predator_prey_ratio` (default `0`)
- `specificity` (default `0`)
- `protection` (default `0`)

Traits may be numeric constants or `AbstractParamDef` instances.

Missing traits default to inactive values; explicit `nothing` is treated the same
as missing (inactive). This matches Agate's fixed missing/nothing semantics.
"""
function build_palatability_matrix(
    ::Type{FT},
    pft_data,
    diameters::AbstractVector{FT};
    palatability_fn=allometric_palatability_unimodal_protection,
) where {FT<:AbstractFloat}
    n = length(diameters)

    can_eat = Vector{Bool}(undef, n)
    can_be_eaten = Vector{Bool}(undef, n)
    optimum = Vector{FT}(undef, n)
    spec = Vector{FT}(undef, n)
    prot = Vector{FT}(undef, n)

    @inbounds for i in 1:n
        pd = pft_data[i]
        can_eat[i] = Bool(_trait_get(pd, :can_eat, false))
        can_be_eaten[i] = Bool(_trait_get(pd, :can_be_eaten, true))
        optimum[i] = _resolve_trait(FT, _trait_get(pd, :optimum_predator_prey_ratio, zero(FT)), diameters[i])
        spec[i] = _resolve_trait(FT, _trait_get(pd, :specificity, zero(FT)), diameters[i])
        prot[i] = _resolve_trait(FT, _trait_get(pd, :protection, zero(FT)), diameters[i])
    end

    return palatability_matrix_allometric(
        FT,
        diameters;
        can_eat=can_eat,
        can_be_eaten=can_be_eaten,
        optimum_predator_prey_ratio=optimum,
        specificity=spec,
        protection=prot,
        palatability_fn=palatability_fn,
    )
end

"""
    build_assimilation_matrix(FT, pft_data, diameters)

Build the default assimilation-efficiency matrix `β[pred, prey]`.

Traits used by the default rule:

- `can_eat::Bool` (default `false`)
- `can_be_eaten::Bool` (default `true`)
- `assimilation_efficiency` (default `0`)
"""
function build_assimilation_matrix(
    ::Type{FT},
    pft_data,
    diameters::AbstractVector{FT};
) where {FT<:AbstractFloat}
    n = length(diameters)

    can_eat = Vector{Bool}(undef, n)
    can_be_eaten = Vector{Bool}(undef, n)
    assim = Vector{FT}(undef, n)

    @inbounds for i in 1:n
        pd = pft_data[i]
        can_eat[i] = Bool(_trait_get(pd, :can_eat, false))
        can_be_eaten[i] = Bool(_trait_get(pd, :can_be_eaten, true))
        assim[i] = _resolve_trait(FT, _trait_get(pd, :assimilation_efficiency, zero(FT)), diameters[i])
    end

    return assimilation_efficiency_matrix_binary(
        FT;
        can_eat=can_eat,
        can_be_eaten=can_be_eaten,
        assimilation_efficiency=assim,
    )
end

end # module Allometry
