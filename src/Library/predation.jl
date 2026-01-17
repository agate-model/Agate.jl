"""Predation and grazing functors."""

module Predation

export HollingTypeII,
    IdealizedPredationLoss,
    IdealizedPredationGain,
    IdealizedPredationAssimilationLoss,
    PreferentialPredationLoss,
    PreferentialPredationGain,
    PreferentialPredationAssimilationLoss,
    AssimilationPreyParameters,
    AssimilationPredatorParameters,
    EmergentAssimilationEfficiencyBinary

"""
    HollingTypeII(K)

Holling type-II saturation factor `P ↦ P / (K + P)`.
"""
struct HollingTypeII{T}
    K::T
end

@inline (h::HollingTypeII)(P) = P / (h.K + P)

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
    loss = PreferentialPredationLoss(f.maximum_grazing_rate, f.half_saturation, f.palatability)(P, Z)
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
    loss = PreferentialPredationLoss(f.maximum_grazing_rate, f.half_saturation, f.palatability)(P, Z)
    return (one(f.assimilation_efficiency) - f.assimilation_efficiency) * loss
end

"""Parameters describing prey traits relevant to assimilation."""
Base.@kwdef struct AssimilationPreyParameters
    can_be_eaten::Bool = false
end

"""Parameters describing predator traits relevant to assimilation."""
Base.@kwdef struct AssimilationPredatorParameters{T}
    can_eat::Bool = false
    assimilation_efficiency::T = 0.0
end

"""
    EmergentAssimilationEfficiencyBinary()

Binary emergent assimilation efficiency controlled by predator/prey traits.
"""
struct EmergentAssimilationEfficiencyBinary end

@inline function (f::EmergentAssimilationEfficiencyBinary)(
    prey::AssimilationPreyParameters,
    predator::AssimilationPredatorParameters,
)
    return ifelse(predator.can_eat && prey.can_be_eaten, predator.assimilation_efficiency, zero(predator.assimilation_efficiency))
end

end # module
