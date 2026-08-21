"""Predation and grazing functors."""

module Predation

export holling_type_ii, idealized_predation_loss, preferential_predation_loss

"""
    HollingTypeII(K)

Holling (1959) type-II functional response factor.

!!! formulation
    ``R`` / (``K`` + ``R``)

    where:
    - ``R`` = prey concentration
    - ``K`` = prey half-saturation constant

The formulation is characterized by decelerating predation as prey concentrations increase.
"""
struct HollingTypeII{T}
    K::T
end

@inline function (h::HollingTypeII)(P)
    K = h.K
    if K == zero(K) && P == zero(P)
        return zero(P)
    end
    return P / (K + P)
end

"""
    PreferentialPredationLoss(maximum_grazing_rate, half_saturation, palatability)

Preferential predation loss from prey `P` to predator `Z`.

!!! formulation
    gₘₐₓ * η * γᴾᴿᴱᴰ * Z

    where:
    - gₘₐₓ = maximum grazing rate
    - η = palatability
    - γᴾᴿᴱᴰ = `HollingTypeII(K)(P)`
    - P = prey concentration
    - K = prey half-saturation
    - Z = predator concentration
"""
struct PreferentialPredationLoss{T1,T2,T3}
    maximum_grazing_rate::T1
    half_saturation::T2
    palatability::T3
end

@inline function (f::PreferentialPredationLoss)(P, Z)
    return f.maximum_grazing_rate * f.palatability * HollingTypeII(f.half_saturation)(P) * Z
end

# -----------------------------------------------------------------------------
# Explicit function aliases (preferred developer UX).
# -----------------------------------------------------------------------------

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
@inline holling_type_ii(P, K) = HollingTypeII(K)(P)

"""
    preferential_predation_loss(P, Z, maximum_grazing_rate, half_saturation, palatability)

Preferential predation loss from prey `P` to predator `Z`.

!!! formulation
    gₘₐₓ * η * (P / (K + P)) * Z

# Arguments
- `P`: prey concentration
- `Z`: predator concentration
- `maximum_grazing_rate`: maximum grazing rate gₘₐₓ
- `half_saturation`: prey half-saturation K
- `palatability`: palatability η
"""
@inline preferential_predation_loss(P, Z, maximum_grazing_rate, half_saturation, palatability) = PreferentialPredationLoss(
    maximum_grazing_rate, half_saturation, palatability
)(
    P, Z
)

end # module
