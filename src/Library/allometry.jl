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
    allometric_scaling_power(a::Number, b::Number, diameter::Number)

Allometric scaling function using the power law for cell volume.

!!! formulation
    ``a````V````ᵇ``

    where:
    - ``V`` = (4 / 3) * π * (``d`` / 2)³
    - ``a`` = scale
    - ``b`` = exponent
    - ``d`` = cell equivalent spherical diameter (ESD)

# Arguments
- `a::Number`: scale parameter.
- `b::Number`: exponent parameter.
- `diameter::Number`: cell equivalent spherical diameter (ESD), in the same length units used throughout the model.

# Returns
- Scaled value `a * V^b` as `FT`, where `V` is the spherical volume computed from `diameter`.
"""
@inline function allometric_scaling_power(a::FT, b::FT, diameter::FT) where {FT<:AbstractFloat}
    r = diameter / FT(2)
    volume = (FT(4) / FT(3)) * FT(π) * r^FT(3)
    return a * volume^b
end

"""
    allometric_palatability_unimodal(prey::PalatabilityPreyParameters, predator::PalatabilityPredatorParameters)

Calculates the unimodal allometric palatability of prey based on predator-prey diameters.

!!! formulation
    0 if ``e_{pred}`` = 0

    1 / (1 + (``d_{ratio}``- ``d_{opt}``)²)``^σ``  otherwise

    where:
    - ``e_{pred}`` = binary ability of predator to eat prey
    - ``d_{ratio}`` = ratio between predator and prey diameters
    - ``d_{opt}`` = optimum ratio between predator and prey diameter
    - σ = how sharply the palatability decreases away from the optimal ratio.

!!! info
    This formulation differs from the operational MITgcm-DARWIN model as it is is structurally different and diameters are used instead of volumes.
    However, both formulations result in a unimodal response where the width and optima are modulated by the optimum-predator-prey ratio and the specificity.

# Arguments
- `prey::PalatabilityPreyParameters{FT}`: prey parameters.
  - `prey.diameter::FT`: prey equivalent spherical diameter (ESD).
  - `prey.protection::FT`: prey protection (unused in this function; see `allometric_palatability_unimodal_protection`).
- `predator::PalatabilityPredatorParameters{FT}`: predator parameters.
  - `predator.can_eat::Bool`: whether the predator can consume this prey type.
  - `predator.diameter::FT`: predator equivalent spherical diameter (ESD).
  - `predator.optimum_predator_prey_ratio::FT`: optimal predator:prey diameter ratio (`d_opt`).
  - `predator.specificity::FT`: unimodal sharpness parameter (σ).

# Returns
- `FT`: palatability in `[0, 1]` (returns `0` when `predator.can_eat` is `false`).
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
    allometric_palatability_unimodal_protection(prey::PalatabilityPreyParameters, predator::PalatabilityPredatorParameters)

Calculates the unimodal allometric palatability of prey, accounting for additional prey protection mechanisms.

!!! formulation
    0 if ``e_{pred}`` = 0

    (1 - η) / (1 + (``d_{ratio}``- ``d_{opt}``)^2)``^σ``   otherwise

    where:
    - ``e_{pred}`` = binary ability of predator to eat prey
    - η = prey-protection
    - ``d_{ratio}`` = ratio between predator and prey diameters
    - ``d_{opt}`` = optimum ratio between predator and prey diameter
    - σ = how sharply the palatability decreases away from the optimal ratio.

# Arguments
- `prey::PalatabilityPreyParameters{FT}`: prey parameters.
  - `prey.diameter::FT`: prey equivalent spherical diameter (ESD).
  - `prey.protection::FT`: prey protection factor η, typically in `[0, 1]`, reducing palatability as `(1 - η)`.
- `predator::PalatabilityPredatorParameters{FT}`: predator parameters.
  - `predator.can_eat::Bool`: whether the predator can consume this prey type.
  - `predator.diameter::FT`: predator equivalent spherical diameter (ESD).
  - `predator.optimum_predator_prey_ratio::FT`: optimal predator:prey diameter ratio (`d_opt`).
  - `predator.specificity::FT`: unimodal sharpness parameter (σ).

# Returns
- `FT`: palatability in `[0, 1]` (returns `0` when `predator.can_eat` is `false`).
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
