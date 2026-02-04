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

Holling type-II saturation factor `P ↦ P / (K + P)`.
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

Idealized predation loss from a prey `P` to a predator `Z` using a squared Holling term.
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

Assimilated predation gain to a predator `Z` from a prey `P`.
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

Unassimilated fraction of predation loss.
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

Preferential predation loss `gmax * palatability * HollingTypeII(half_saturation)(P) * Z`.
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

Unassimilated fraction of preferential predation loss.
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

"""Holling type-II saturation factor `P/(K+P)`."""
@inline holling_type_ii(P, K) = HollingTypeII(K)(P)

"""Idealized predation loss from prey `P` to predator `Z` (squared Holling)."""
@inline idealized_predation_loss(P, Z, maximum_grazing_rate, half_saturation) =
    IdealizedPredationLoss(maximum_grazing_rate, half_saturation)(P, Z)

"""Assimilated idealized predation gain to predator `Z` from prey `P`."""
@inline idealized_predation_gain(
    P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation
) = IdealizedPredationGain(assimilation_efficiency, maximum_grazing_rate, half_saturation)(
    P, Z
)

"""Unassimilated fraction of idealized predation loss."""
@inline idealized_predation_unassimilated_loss(
    P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation
) = IdealizedPredationAssimilationLoss(
    assimilation_efficiency, maximum_grazing_rate, half_saturation
)(
    P, Z
)

"""Preferential predation loss from prey `P` to predator `Z`."""
@inline preferential_predation_loss(
    P, Z, maximum_grazing_rate, half_saturation, palatability
) = PreferentialPredationLoss(maximum_grazing_rate, half_saturation, palatability)(P, Z)

"""Assimilated preferential predation gain to predator `Z` from prey `P`."""
@inline preferential_predation_gain(
    P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability
) = PreferentialPredationGain(
    assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability
)(
    P, Z
)

"""Unassimilated fraction of preferential predation loss."""
@inline preferential_predation_unassimilated_loss(
    P, Z, assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability
) = PreferentialPredationAssimilationLoss(
    assimilation_efficiency, maximum_grazing_rate, half_saturation, palatability
)(
    P, Z
)

end # module
