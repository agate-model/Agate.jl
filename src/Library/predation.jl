"""Predation and grazing formulations."""

module Predation

export holling_type_ii, idealized_predation_loss, preferential_predation_loss

"""
    idealized_predation_loss(P, Z, maximum_grazing_rate, half_saturation)

Idealized NPZD grazing loss with a squared Holling response.
"""
@inline function idealized_predation_loss(P, Z, maximum_grazing_rate, half_saturation)
    P_squared = P * P
    K_squared = half_saturation * half_saturation
    K_squared == zero(K_squared) && P_squared == zero(P_squared) && return zero(P * Z)
    return maximum_grazing_rate * (P_squared / (K_squared + P_squared)) * Z
end

"""
    holling_type_ii(P, K)

Holling (1959) type-II functional response.

!!! formulation
    ``P`` / (``K`` + ``P``)

# Arguments
- `P`: prey concentration
- `K`: prey half-saturation (prey density at which predation is half its maximum)
"""
@inline function holling_type_ii(P, K)
    K == zero(K) && P == zero(P) && return zero(P)
    return P / (K + P)
end

"""
    preferential_predation_loss(P, Z, maximum_grazing_rate, half_saturation, palatability)

Preferential predation loss from prey `P` to predator `Z`.

!!! formulation
    gmax * eta * (P / (K + P)) * Z

# Arguments
- `P`: prey concentration
- `Z`: predator concentration
- `maximum_grazing_rate`: maximum grazing rate
- `half_saturation`: prey half-saturation
- `palatability`: palatability
"""
@inline function preferential_predation_loss(
    P, Z, maximum_grazing_rate, half_saturation, palatability
)
    return maximum_grazing_rate * palatability * holling_type_ii(P, half_saturation) * Z
end

end # module
