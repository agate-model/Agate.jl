"""Predation and grazing kernels."""
module Predation

export holling_type_ii, linear_predation_loss, proportional_predation_loss, switching_predation_loss

"""
    holling_type_ii(P, K)

Return the Holling (1959) type-II functional response ``P / (K + P)``.
The indeterminate `P == K == 0` case returns zero.
"""
@inline function holling_type_ii(P, K)
    K == zero(K) && P == zero(P) && return zero(P)
    return P / (K + P)
end

"""
    linear_predation_loss(inventory, consumer, grazing_rate, palatability)

Return mass-action prey loss ``g p R Z`` for grazing coefficient `g`, palatability `p`,
prey inventory `R`, and consumer biomass `Z`.
"""
@inline linear_predation_loss(inventory, consumer, grazing_rate, palatability) =
    grazing_rate * palatability * inventory * consumer

"""
    proportional_predation_loss(
        inventory, consumer, maximum_grazing_rate, half_saturation,
        palatability, total_palatable_biomass
    )

Return proportional-allocation prey loss for shared-capacity grazing.
`total_palatable_biomass` is the consumer-level sum ``sum(p_j R_j)`` supplied as a scalar
runtime-IR reduction.
"""
@inline function proportional_predation_loss(
    inventory,
    consumer,
    maximum_grazing_rate,
    half_saturation,
    palatability,
    total_palatable_biomass,
)
    half_saturation == zero(half_saturation) &&
        total_palatable_biomass == zero(total_palatable_biomass) &&
        return zero(maximum_grazing_rate * inventory * consumer)
    return maximum_grazing_rate * palatability * inventory /
           (half_saturation + total_palatable_biomass) * consumer
end

"""
    switching_predation_loss(
        inventory, reference_inventory, consumer, maximum_grazing_rate,
        half_saturation, palatability, total_palatable_biomass,
        switching_weight_sum, switching_exponent
    )

Return switching prey loss for shared-capacity grazing. The two consumer-level reductions are
provided as scalar runtime-IR nodes so edge expressions stay compact as food-web size increases.
"""
@inline function switching_predation_loss(
    inventory,
    reference_inventory,
    consumer,
    maximum_grazing_rate,
    half_saturation,
    palatability,
    total_palatable_biomass,
    switching_weight_sum,
    switching_exponent,
)
    reference_inventory == zero(reference_inventory) &&
        return zero(maximum_grazing_rate * inventory * consumer)
    saturation = holling_type_ii(total_palatable_biomass, half_saturation)
    switching_weight_sum == zero(switching_weight_sum) &&
        return zero(maximum_grazing_rate * inventory * consumer)
    allocation = (palatability * reference_inventory)^switching_exponent / switching_weight_sum
    return maximum_grazing_rate * consumer * saturation * allocation *
           inventory / reference_inventory
end


end # module
