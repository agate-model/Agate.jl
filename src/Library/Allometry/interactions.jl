# -----------------------------------------------------------------------------
# Palatability + assimilation matrix utilities
# -----------------------------------------------------------------------------

"""Prey traits used for allometric palatability."""
struct PalatabilityPreyParameters{FT<:AbstractFloat}
    diameter::FT
    protection::FT
end

"""Predator traits used for allometric palatability."""
struct PalatabilityPredatorParameters{FT<:AbstractFloat}
    diameter::FT
    optimum_predator_prey_ratio::FT
    specificity::FT
end

"""
    allometric_scaling_power(a, b, diameter)

Allometric scaling function using a power law on spherical cell volume.

!!! formulation
    ``a````V````ᵇ``

    where:
    - ``V`` = (4 / 3) * π * (``d`` / 2)³
    - ``a`` = scale
    - ``b`` = exponent
    - ``d`` = cell equivalent spherical diameter (ESD)

# Arguments
- `a`: scale parameter
- `b`: exponent parameter
- `diameter`: cell equivalent spherical diameter (ESD)

# Returns
`a * V^b` where `V` is the spherical volume computed from `diameter`.
"""
@inline function allometric_scaling_power(
    a::FT, b::FT, diameter::FT
) where {FT<:AbstractFloat}
    r = diameter / FT(2)
    volume = (FT(4) / FT(3)) * FT(π) * r^FT(3)
    return a * volume^b
end

"""
    allometric_palatability_unimodal(prey, predator)

Unimodal allometric palatability based on predator–prey diameters.

!!! formulation
    1 / (1 + (``ratio`` - ``opt``)²)``^σ``

    where:
    - ``ratio`` = predator / prey diameter ratio
    - ``opt`` = optimum predator:prey diameter ratio
    - σ = unimodal sharpness parameter (specificity)

# Arguments
- `prey`: `PalatabilityPreyParameters(diameter, protection)`
- `predator`: `PalatabilityPredatorParameters(diameter, optimum_predator_prey_ratio, specificity)`

# Returns
A palatability value in `[0, 1]`.
"""
@inline function allometric_palatability_unimodal(
    prey::PalatabilityPreyParameters{FT}, predator::PalatabilityPredatorParameters{FT}
) where {FT<:AbstractFloat}
    ratio = predator.diameter / prey.diameter
    width = one(FT) + (ratio - predator.optimum_predator_prey_ratio)^FT(2)
    return one(FT) / width^predator.specificity
end

"""
    allometric_palatability_unimodal_protection(prey, predator)

Unimodal allometric palatability with prey protection.

Protection η reduces palatability as `(1 - η)` (so η=0 means no protection).

!!! formulation
    (1 - η) / (1 + (``ratio`` - ``opt``)²)``^σ``

    where η is `prey.protection`.
"""
function allometric_palatability_unimodal_protection(
    prey::PalatabilityPreyParameters{FT}, predator::PalatabilityPredatorParameters{FT}
) where {FT<:AbstractFloat}
    base = allometric_palatability_unimodal(prey, predator)
    return base * (one(FT) - prey.protection)
end

"""
    palatability_matrix_allometric_axes(FT, diameters; optimum_predator_prey_ratio, specificity, protection,
                                        consumer_indices, prey_indices; palatability_fn=allometric_palatability_unimodal_protection)

Build a role-aware palatability matrix `M[consumer, prey]` using allometric traits.

Only the consumer-by-prey block specified by `consumer_indices` and `prey_indices` is
computed, producing a rectangular (n_consumer × n_prey) matrix.
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

"""
    assimilation_efficiency_matrix_binary_axes(FT; assimilation_efficiency, consumer_indices, prey_indices)

Build a role-aware assimilation-efficiency matrix `β[consumer, prey]`.

Only the consumer-by-prey block specified by `consumer_indices` and `prey_indices`
is computed.

The returned matrix assigns each consumer's scalar `assimilation_efficiency` across all
prey indices.
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
