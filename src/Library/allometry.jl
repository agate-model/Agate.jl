module Allometry

"""Utilities for size-dependent traits and interaction matrices."""

export AbstractParamDef, ConstantParam, AllometricParam
export PowerLaw, resolve_param, cast_paramdef

export PalatabilityPreyParameters, PalatabilityPredatorParameters
export allometric_scaling_power
export allometric_palatability_unimodal, allometric_palatability_unimodal_protection
export palatability_matrix_allometric_axes, assimilation_efficiency_matrix_binary_axes
# -----------------------------------------------------------------------------
# Explicit parameter definitions
# -----------------------------------------------------------------------------

"""Abstract supertype for explicit parameter definitions stored in `PFTSpecification`."""
abstract type AbstractParamDef end

"""A parameter that is constant across sizes."""
struct ConstantParam{T} <: AbstractParamDef
    value::T
end

# NOTE: the default outer constructor `ConstantParam(x)` already exists and infers `T`.
# We avoid redefining it to prevent method overwrite warnings during precompilation.

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
- `prefactor`
- `exponent`

(Backward-compat coefficient aliases were intentionally removed.)
"""
@inline function (m::PowerLaw)(coeffs::NamedTuple, diameter)
    hasproperty(coeffs, :prefactor) ||
        throw(ArgumentError("PowerLaw requires coefficient `prefactor`"))
    hasproperty(coeffs, :exponent) ||
        throw(ArgumentError("PowerLaw requires coefficient `exponent`"))

    a = getproperty(coeffs, :prefactor)
    b = getproperty(coeffs, :exponent)

    # By construction we keep coefficients and diameter the same floating type (FT),
    # so this call never mixes Float32/Float64 (important for GPU use).
    return allometric_scaling_power(a, b, diameter)
end

"""Resolve a parameter definition at a given diameter.

This is the one function the constructor uses when building runtime parameter vectors.
"""
@inline resolve_param(::Type{FT}, x, diameter) where {FT<:AbstractFloat} = FT(x)

@inline resolve_param(::Type{FT}, x::Bool, diameter) where {FT<:AbstractFloat} = x

@inline resolve_param(::Type{FT}, p::ConstantParam, diameter) where {FT<:AbstractFloat} =
    FT(p.value)

@inline function resolve_param(
    ::Type{FT}, p::AllometricParam, diameter
) where {FT<:AbstractFloat}
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
    diameter::FT
    optimum_predator_prey_ratio::FT
    specificity::FT
end

"""Power-law scaling function using spherical volume."""
@inline function allometric_scaling_power(
    a::FT, b::FT, diameter::FT
) where {FT<:AbstractFloat}
    r = diameter / FT(2)
    volume = (FT(4) / FT(3)) * FT(π) * r^FT(3)
    return a * volume^b
end

"""Unimodal palatability (no protection)."""
@inline function allometric_palatability_unimodal(
    prey::PalatabilityPreyParameters{FT}, predator::PalatabilityPredatorParameters{FT}
) where {FT<:AbstractFloat}
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

"""Build a role-aware allometric palatability matrix `M[consumer, prey]`.

Only the consumer-by-prey block specified by `consumer_indices` and `prey_indices`
is computed, producing a rectangular (n_consumer × n_prey) matrix.
"""
function palatability_matrix_allometric_axes(
    ::Type{FT},
    diameters::AbstractVector{FT};
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
        predator = PalatabilityPredatorParameters{FT}(
            diameters[pred], optimum_predator_prey_ratio[pred], specificity[pred]
        )
        for (jj, prey) in pairs(prey_indices)
            prey_params = PalatabilityPreyParameters{FT}(diameters[prey], protection[prey])
            M[ii, jj] = palatability_fn(prey_params, predator)
        end
    end

    return M
end

"""Build a role-aware binary assimilation-efficiency matrix `β[consumer, prey]`.

Only the consumer-by-prey block specified by `consumer_indices` and `prey_indices`
is computed.
"""
function assimilation_efficiency_matrix_binary_axes(
    ::Type{FT};
    assimilation_efficiency::AbstractVector{FT},
    consumer_indices::AbstractVector{<:Integer},
    prey_indices::AbstractVector{<:Integer},
) where {FT<:AbstractFloat}
    nr = length(consumer_indices)
    nc = length(prey_indices)
    M = zeros(FT, nr, nc)

    @inbounds for (ii, pred) in pairs(consumer_indices)
        β = assimilation_efficiency[pred]
        for jj in 1:nc
            M[ii, jj] = β
        end
    end

    return M
end

end # module Allometry
