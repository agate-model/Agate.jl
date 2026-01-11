"""NiPiZD defaults expressed as a model-agnostic factory.

This file defines `NiPiZDFactory` and the default inputs used by
`Agate.Models.construct(factory; ...)`.
"""

using Oceananigans.Units

using Agate.Utils: AbstractBGCFactory, PFTParameters, BiogeochemistrySpecification
using Agate.Utils: DiameterRangeSpecification
using Agate.Utils: InteractionDynamics
using Agate.Utils: pft_get
import Agate.Utils: default_interactions

import Agate.Models: default_plankton_dynamics, default_plankton_args, default_biogeochem_dynamics, default_biogeochem_args
using .Tracers:
    nutrient_default,
    detritus_default,
    phytoplankton_default,
    zooplankton_default

using Agate.Library.Allometry:
    PalatabilityPreyParameters,
    PalatabilityPredatorParameters,
    allometric_palatability_unimodal_protection

using Agate.Library.Predation:
    AssimilationPreyParameters,
    AssimilationPredatorParameters,
    assimilation_efficiency_emergent_binary

"""Factory for the size-structured NiPiZD model."""
struct NiPiZDFactory <: AbstractBGCFactory end

"""Default plankton dynamics for NiPiZD.

Returns a `NamedTuple` mapping group prefix => tracer dynamics builder.
"""
default_plankton_dynamics(::NiPiZDFactory) = (
    Z = zooplankton_default,
    P = phytoplankton_default,
)

"""Default plankton arguments for NiPiZD.

Returns a `NamedTuple` mapping group prefix => group specification.

Ordering is significant; the default keeps the historical `Z`-then-`P` ordering.
"""
function default_plankton_args(::NiPiZDFactory, ::Type{FT}) where {FT<:AbstractFloat}
    phyto_pft = PFTParameters(
        maximum_growth_rate_a = FT(2 / day),
        maximum_growth_rate_b = FT(-0.15),
        nutrient_half_saturation_a = FT(0.17),
        nutrient_half_saturation_b = FT(0.27),
        linear_mortality = FT(8e-7 / second),
        alpha = FT(0.1953 / day),
        photosynthetic_slope = zero(FT),
        chlorophyll_to_carbon_ratio = zero(FT),
        can_eat = false,
        can_be_eaten = true,
        optimum_predator_prey_ratio = zero(FT),
        protection = zero(FT),
        specificity = zero(FT),
        assimilation_efficiency = zero(FT),
    )

    zoo_pft = PFTParameters(
        maximum_predation_rate_a = FT(30.84 / day),
        maximum_predation_rate_b = FT(-0.16),
        linear_mortality = FT(8e-7 / second),
        holling_half_saturation = FT(5.0),
        quadratic_mortality = FT(1e-6 / second),
        can_eat = true,
        can_be_eaten = false,
        optimum_predator_prey_ratio = FT(10),
        protection = one(FT),
        specificity = FT(0.3),
        assimilation_efficiency = FT(0.32),
    )

    return (
        Z = (; n = 2, diameters = DiameterRangeSpecification(20, 100, :linear_splitting), pft = zoo_pft),
        P = (; n = 2, diameters = DiameterRangeSpecification(2, 10, :log_splitting), pft = phyto_pft),
    )
end

"""Default non-plankton tracer dynamics for NiPiZD."""
default_biogeochem_dynamics(::NiPiZDFactory) = (
    N = nutrient_default,
    D = detritus_default,
)

"""Default biogeochemical specification for NiPiZD."""
function default_biogeochem_args(::NiPiZDFactory, ::Type{FT}) where {FT<:AbstractFloat}
    return BiogeochemistrySpecification(
        detritus_remineralization = FT(0.1213 / day),
        mortality_export_fraction = FT(0.2),
    )
end

"""Default interaction builder for NiPiZD.

Computes palatability and assimilation-efficiency matrices using allometry.
"""
function default_interactions(::NiPiZDFactory)
    return InteractionDynamics(_nipizd_default_interactions)
end

function _nipizd_default_interactions(ctx, args=NamedTuple())
    FT = ctx.FT
    n = ctx.n_total
    diameters = ctx.diameters
    pfts = ctx.pfts

    pal = zeros(FT, n, n)
    assim = zeros(FT, n, n)

    @inbounds for pred in 1:n
        pred_pft = pfts[pred]
        predator = PalatabilityPredatorParameters{FT}(
            Bool(pft_get(pred_pft, :can_eat, false)),
            diameters[pred],
            FT(pft_get(pred_pft, :optimum_predator_prey_ratio, zero(FT))),
            FT(pft_get(pred_pft, :specificity, zero(FT))),
        )

        predator_assim = AssimilationPredatorParameters{FT}(
            Bool(pft_get(pred_pft, :can_eat, false)),
            FT(pft_get(pred_pft, :assimilation_efficiency, zero(FT))),
        )

        for prey in 1:n
            prey_pft = pfts[prey]

            prey_params = PalatabilityPreyParameters{FT}(
                diameters[prey],
                FT(pft_get(prey_pft, :protection, zero(FT))),
            )
            pal[pred, prey] = allometric_palatability_unimodal_protection(prey_params, predator)

            prey_assim = AssimilationPreyParameters(Bool(pft_get(prey_pft, :can_be_eaten, true)))
            assim[pred, prey] = assimilation_efficiency_emergent_binary(prey_assim, predator_assim)
        end
    end

    return (; palatability_matrix = pal, assimilation_efficiency_matrix = assim)
end
