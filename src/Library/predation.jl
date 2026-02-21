"""Predation and grazing functors."""

module Predation

export holling_type_ii,
    idealized_predation_loss,
    idealized_predation_gain,
    idealized_predation_unassimilated_loss,
    preferential_predation_loss,
    preferential_predation_gain,
    preferential_predation_unassimilated_loss

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
    IdealizedPredationLoss(maximum_grazing_rate, half_saturation)

Idealized loss rate of prey `P` to a predator `Z` using a squared Holling term.

!!! formulation
    gₘₐₓ * γᴾᴿᴱᴰ * Z

    where:
    - gₘₐₓ = maximum grazing rate of the predator
    - γᴾᴿᴱᴰ = `P² / (K² + P²)`
    - P = prey concentration
    - K = prey half-saturation
    - Z = predator concentration
"""
struct IdealizedPredationLoss{T}
    maximum_grazing_rate::T
    half_saturation::T
end

@inline function (f::IdealizedPredationLoss)(P, Z)
    gmax = f.maximum_grazing_rate
    K = f.half_saturation
    P2 = P * P
    if K == zero(K) && P2 == zero(P2)
        return zero(P * Z)
    end
    return gmax * (P2 / (K * K + P2)) * Z
end

"""
    IdealizedPredationGain(assimilation_efficiency, maximum_grazing_rate, half_saturation)

Assimilated gain rate to predator `Z` from prey `P`.

!!! formulation
    β * g

    where:
    - β = assimilation efficiency
    - g = `IdealizedPredationLoss(gₘₐₓ, K)(P, Z)`
"""
struct IdealizedPredationGain{T}
    assimilation_efficiency::T
    maximum_grazing_rate::T
    half_saturation::T
end

@inline function (f::IdealizedPredationGain)(P, Z)
    loss = IdealizedPredationLoss(f.maximum_grazing_rate, f.half_saturation)(P, Z)
    return f.assimilation_efficiency * loss
end

"""
    IdealizedPredationAssimilationLoss(assimilation_efficiency, maximum_grazing_rate, half_saturation)

Unassimilated fraction of idealized predation loss ("sloppy feeding").

!!! formulation
    (1 - β) * g

    where:
    - β = assimilation efficiency
    - g = `IdealizedPredationLoss(gₘₐₓ, K)(P, Z)`

!!! note
    This represents transfer of biomass from prey to the environment rather than to the predator.
"""
struct IdealizedPredationAssimilationLoss{T}
    assimilation_efficiency::T
    maximum_grazing_rate::T
    half_saturation::T
end

@inline function (f::IdealizedPredationAssimilationLoss)(P, Z)
    loss = IdealizedPredationLoss(f.maximum_grazing_rate, f.half_saturation)(P, Z)
    return (one(f.assimilation_efficiency) - f.assimilation_efficiency) * loss
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
struct PreferentialPredationLoss{T}
    maximum_grazing_rate::T
    half_saturation::T
    palatability::T
end

@inline function (f::PreferentialPredationLoss)(P, Z)
    return f.maximum_grazing_rate * f.palatability * HollingTypeII(f.half_saturation)(P) * Z
end

"""
    PreferentialPredationGain(assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability)

Assimilated preferential predation gain.

!!! formulation
    β * g

    where:
    - β = assimilation efficiency
    - g = `PreferentialPredationLoss(gₘₐₓ, K, η)(P, Z)`
"""
struct PreferentialPredationGain{T}
    assimilation_efficiency::T
    maximum_grazing_rate::T
    half_saturation::T
    palatability::T
end

@inline function (f::PreferentialPredationGain)(P, Z)
    loss = PreferentialPredationLoss(
        f.maximum_grazing_rate, f.half_saturation, f.palatability
    )(
        P, Z
    )
    return f.assimilation_efficiency * loss
end

"""
    PreferentialPredationAssimilationLoss(assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability)

Unassimilated fraction of preferential predation loss ("sloppy feeding").

!!! formulation
    (1 - β) * g

    where:
    - β = assimilation efficiency
    - g = `PreferentialPredationLoss(gₘₐₓ, K, η)(P, Z)`

!!! info
    This represents transfer of biomass from prey to the environment rather than to the predator.
"""
struct PreferentialPredationAssimilationLoss{T}
    assimilation_efficiency::T
    maximum_grazing_rate::T
    half_saturation::T
    palatability::T
end

@inline function (f::PreferentialPredationAssimilationLoss)(P, Z)
    loss = PreferentialPredationLoss(
        f.maximum_grazing_rate, f.half_saturation, f.palatability
    )(
        P, Z
    )
    return (one(f.assimilation_efficiency) - f.assimilation_efficiency) * loss
end

# -----------------------------------------------------------------------------
# Explicit function aliases (preferred developer UX).
# -----------------------------------------------------------------------------

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
    idealized_predation_loss(P, Z, maximum_grazing_rate, half_saturation)

Loss rate of prey `P` to predator `Z` using a squared Holling term.

!!! formulation
    gₘₐₓ * (P² / (K² + P²)) * Z

# Arguments
- `P`: prey concentration
- `Z`: predator concentration
- `maximum_grazing_rate`: maximum grazing rate gₘₐₓ
- `half_saturation`: prey half-saturation K
"""
@inline idealized_predation_loss(P, Z, maximum_grazing_rate, half_saturation) =
    IdealizedPredationLoss(maximum_grazing_rate, half_saturation)(P, Z)

"""
    idealized_predation_gain(P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation)

Assimilated gain rate to predator `Z` feeding on prey `P`.

!!! formulation
    β * g

    where `g = idealized_predation_loss(P, Z, gₘₐₓ, K)`.

# Arguments
- `P`: prey concentration
- `Z`: predator concentration
- `assimilation_efficiency`: assimilation efficiency β
- `maximum_grazing_rate`: maximum grazing rate gₘₐₓ
- `half_saturation`: prey half-saturation K
"""
@inline idealized_predation_gain(
    P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation
) = IdealizedPredationGain(assimilation_efficiency, maximum_grazing_rate, half_saturation)(
    P, Z
)

"""
    idealized_predation_unassimilated_loss(P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation)

Unassimilated fraction of idealized predation loss ("sloppy feeding").

!!! formulation
    (1 - β) * g

    where `g = idealized_predation_loss(P, Z, gₘₐₓ, K)`.

# Arguments
- `P`: prey concentration
- `Z`: predator concentration
- `assimilation_efficiency`: assimilation efficiency β
- `maximum_grazing_rate`: maximum grazing rate gₘₐₓ
- `half_saturation`: prey half-saturation K
"""
@inline idealized_predation_unassimilated_loss(
    P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation
) = IdealizedPredationAssimilationLoss(
    assimilation_efficiency, maximum_grazing_rate, half_saturation
)(
    P, Z
)

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
@inline preferential_predation_loss(
    P, Z, maximum_grazing_rate, half_saturation, palatability
) = PreferentialPredationLoss(maximum_grazing_rate, half_saturation, palatability)(P, Z)

"""
    preferential_predation_gain(P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability)

Assimilated preferential predation gain.

!!! formulation
    β * g

    where `g = preferential_predation_loss(P, Z, gₘₐₓ, K, η)`.

# Arguments
- `P`: prey concentration
- `Z`: predator concentration
- `assimilation_efficiency`: assimilation efficiency β
- `maximum_grazing_rate`: maximum grazing rate gₘₐₓ
- `half_saturation`: prey half-saturation K
- `palatability`: palatability η
"""
@inline preferential_predation_gain(
    P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability
) = PreferentialPredationGain(
    assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability
)(
    P, Z
)

"""
    preferential_predation_unassimilated_loss(P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability)

Unassimilated fraction of preferential predation loss ("sloppy feeding").

!!! formulation
    (1 - β) * g

    where `g = preferential_predation_loss(P, Z, gₘₐₓ, K, η)`.

# Arguments
- `P`: prey concentration
- `Z`: predator concentration
- `assimilation_efficiency`: assimilation efficiency β
- `maximum_grazing_rate`: maximum grazing rate gₘₐₓ
- `half_saturation`: prey half-saturation K
- `palatability`: palatability η
"""
@inline preferential_predation_unassimilated_loss(
    P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability
) = PreferentialPredationAssimilationLoss(
    assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability
)(
    P, Z
)

end # module
