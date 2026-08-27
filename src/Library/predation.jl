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
    preferential_predation_loss(P, Z, maximum_grazing_rate, half_saturation, palatability)

Return preferential predation loss from prey `P` to predator `Z`:
`maximum_grazing_rate * palatability * holling_type_ii(P, half_saturation) * Z`.
"""
@inline function preferential_predation_loss(
    P, Z, maximum_grazing_rate, half_saturation, palatability
)
    return maximum_grazing_rate * palatability * holling_type_ii(P, half_saturation) * Z
end

end # module
