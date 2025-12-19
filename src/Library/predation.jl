"""Functions related to predator-prey interactions.

This module provides scalar formulations for grazing losses and gains. Functions are
allocation-free and safe to call from GPU kernels.

Agate uses micro-parameter structs for interaction inputs that are passed into
allocation-free scalar functions.
"""

module Predation

export AssimilationPreyParameters
export AssimilationPredatorParameters
export holling_type_2
export predation_loss_idealized
export predation_gain_idealized
export predation_assimilation_loss_idealized
export predation_loss_preferential
export predation_gain_preferential
export predation_assimilation_loss_preferential
export assimilation_efficiency_emergent_binary

"""
    holling_type_2(prey_concentration, prey_half_saturation)

Holling's type II functional response:

```math
R / (k + R)
```
"""
@inline function holling_type_2(prey_concentration, prey_half_saturation)
    return prey_concentration / (prey_half_saturation + prey_concentration)
end

"""
    predation_loss_idealized(P, Z, maximum_grazing_rate, prey_half_saturation)

Loss rate of prey `P` to predator `Z` using an idealized Holling type II response.
"""
@inline function predation_loss_idealized(P, Z, maximum_grazing_rate, prey_half_saturation)
    return maximum_grazing_rate * holling_type_2(P^2, prey_half_saturation^2) * Z
end

"""
    predation_gain_idealized(P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation)

Predator gain rate from feeding on prey with assimilation efficiency.
"""
@inline function predation_gain_idealized(P, Z, assimilation_efficiency, maximum_grazing_rate, k_p)
    return assimilation_efficiency * predation_loss_idealized(P, Z, maximum_grazing_rate, k_p)
end

"""
    predation_assimilation_loss_idealized(P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation)

Loss to the environment due to imperfect assimilation (sloppy feeding).
"""
@inline function predation_assimilation_loss_idealized(P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation)
    return (1 - assimilation_efficiency) * predation_loss_idealized(P, Z, maximum_grazing_rate, prey_half_saturation)
end

"""
    predation_loss_preferential(P, Z, maximum_grazing_rate, prey_half_saturation, palatability)

Loss rate of prey `P` to predator `Z` modulated by palatability.
"""
@inline function predation_loss_preferential(P, Z, maximum_grazing_rate, prey_half_saturation, palatability)
    return maximum_grazing_rate * palatability * holling_type_2(P, prey_half_saturation) * Z
end

"""
    predation_gain_preferential(P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation, palatability)

Predator gain rate from feeding on prey with preferential palatability.
"""
@inline function predation_gain_preferential(P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation, palatability)
    return assimilation_efficiency * predation_loss_preferential(P, Z, maximum_grazing_rate, prey_half_saturation, palatability)
end

"""
    predation_assimilation_loss_preferential(P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation, palatability)

Loss to the environment from preferential predation due to imperfect assimilation.
"""
@inline function predation_assimilation_loss_preferential(P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation, palatability)
    return (1 - assimilation_efficiency) * predation_loss_preferential(P, Z, maximum_grazing_rate, prey_half_saturation, palatability)
end

"""Binary prey data for emergent assimilation efficiency."""
struct AssimilationPreyParameters
    can_be_eaten::Bool
end

"""Predator data for emergent assimilation efficiency."""
struct AssimilationPredatorParameters{FT<:AbstractFloat}
    can_eat::Bool
    assimilation_efficiency::FT
end

"""
    assimilation_efficiency_emergent_binary(prey, predator)

Return the predator assimilation efficiency when the predator can eat and the prey can be eaten.
Returns zero otherwise.
"""
@inline function assimilation_efficiency_emergent_binary(prey::AssimilationPreyParameters, predator::AssimilationPredatorParameters{FT}) where {FT<:AbstractFloat}
    if predator.can_eat && prey.can_be_eaten
        return predator.assimilation_efficiency
    end
    return zero(FT)
end

end # module
