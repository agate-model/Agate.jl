"""Allometric relationships and predator-prey palatability.

This module provides allocation-free scalar functions that are safe to call from GPU
kernels.

Agate uses micro-parameter structs for predator-prey interaction inputs to ensure
interaction functions are type-stable and do not depend on dynamic containers.
"""

module Allometry

export PalatabilityPreyParameters
export PalatabilityPredatorParameters
export allometric_scaling_power
export allometric_palatability_unimodal
export allometric_palatability_unimodal_protection

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

"""
    allometric_scaling_power(a, b, diameter)

Compute an allometrically scaled rate using a power law in cell volume.

The volume is computed assuming an equivalent spherical diameter. The returned value is
``a * V^b``, where ``V = (4/3)π(d/2)^3``.
"""
@inline function allometric_scaling_power(a::FT, b::FT, diameter::FT) where {FT<:AbstractFloat}
    r = diameter / FT(2)
    volume = (FT(4) / FT(3)) * FT(π) * r^FT(3)
    return a * volume^b
end

"""
    allometric_palatability_unimodal(prey, predator)

Unimodal palatability as a function of predator-to-prey diameter ratio.

Returns zero when `predator.can_eat` is false.
"""
@inline function allometric_palatability_unimodal(
    prey::PalatabilityPreyParameters{FT},
    predator::PalatabilityPredatorParameters{FT},
) where {FT<:AbstractFloat}
    predator.can_eat || return zero(FT)
    ratio = predator.diameter / prey.diameter
    width = one(FT) + (ratio - predator.optimum_predator_prey_ratio)^FT(2)
    return one(FT) / width^predator.specificity
end

"""
    allometric_palatability_unimodal_protection(prey, predator)

Unimodal palatability with additional prey protection.
"""
@inline function allometric_palatability_unimodal_protection(
    prey::PalatabilityPreyParameters{FT},
    predator::PalatabilityPredatorParameters{FT},
) where {FT<:AbstractFloat}
    predator.can_eat || return zero(FT)
    ratio = predator.diameter / prey.diameter
    width = one(FT) + (ratio - predator.optimum_predator_prey_ratio)^FT(2)
    return (one(FT) - prey.protection) / width^predator.specificity
end

end # module
