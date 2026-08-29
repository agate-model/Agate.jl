"""Predation and grazing kernels."""
module Predation

export holling_type_ii, preferential_predation_loss

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
    preferential_predation_loss(
        inventory, reference_inventory, consumer, maximum_grazing_rate,
        half_saturation, palatability
    )

Return preferential grazing loss from one prey-state `inventory`, using
`reference_inventory` to determine the shared grazing intensity.
"""
@inline function preferential_predation_loss(
    inventory, reference_inventory, consumer, maximum_grazing_rate, half_saturation, palatability
)
    half_saturation == zero(half_saturation) && reference_inventory == zero(reference_inventory) &&
        return zero(maximum_grazing_rate * inventory * consumer)
    return maximum_grazing_rate * palatability * inventory /
           (half_saturation + reference_inventory) * consumer
end

end # module
