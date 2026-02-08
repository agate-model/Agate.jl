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
