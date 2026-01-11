"""DARWIN defaults expressed as a model-agnostic factory.

This file defines `DarwinFactory` and the default inputs used by
`Agate.Constructor.construct(factory; ...)`.
"""

using Oceananigans.Units

using Agate.Utils: AbstractBGCFactory
using Agate.Utils.Specifications: PFTSpecification, BiogeochemistrySpecification, pft_get
using Agate.Utils: DiameterRangeSpecification
using Agate.Utils: InteractionDynamics, pft_get
import Agate.Utils: default_interactions

import Agate.Models: default_plankton_dynamics, default_plankton_args, default_biogeochem_dynamics, default_biogeochem_args
using .Tracers:
    DIC_geider_light,
    DIN_geider_light,
    PO4_geider_light,
    DOC_default,
    POC_default,
    DON_default,
    PON_default,
    DOP_default,
    POP_default,
    phytoplankton_growth_two_nutrients_geider_light,
    zooplankton_default

using Agate.Library.Allometry:
    PalatabilityPreyParameters,
    PalatabilityPredatorParameters,
    allometric_palatability_unimodal_protection

using Agate.Library.Predation:
    AssimilationPreyParameters,
    AssimilationPredatorParameters,
    assimilation_efficiency_emergent_binary

"""Factory for the simplified DARWIN-like elemental cycling model."""
struct DarwinFactory <: AbstractBGCFactory end

"""Default plankton dynamics for DARWIN."""
default_plankton_dynamics(::DarwinFactory) = (
    Z = zooplankton_default,
    P = phytoplankton_growth_two_nutrients_geider_light,
)

"""Default plankton arguments for DARWIN.

Ordering is significant; the default keeps the historical `Z`-then-`P` ordering.
"""
function default_plankton_args(::DarwinFactory, ::Type{FT}) where {FT<:AbstractFloat}
    phyto_pft = PFTSpecification(
        maximum_growth_rate_a = FT(2 / day),
        maximum_growth_rate_b = FT(-0.15),
        half_saturation_DIN_a = FT(0.17),
        half_saturation_DIN_b = FT(0.27),
        half_saturation_PO4_a = FT(0.17),
        half_saturation_PO4_b = FT(0.27),
        linear_mortality = FT(8e-7 / second),
        photosynthetic_slope = FT(0.46e-5),
        chlorophyll_to_carbon_ratio = FT(0.1),
        can_eat = false,
        can_be_eaten = true,
        optimum_predator_prey_ratio = zero(FT),
        protection = zero(FT),
        specificity = zero(FT),
        assimilation_efficiency = zero(FT),
        alpha = zero(FT),
    )

    zoo_pft = PFTSpecification(
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

"""Default non-plankton tracer dynamics for DARWIN."""
default_biogeochem_dynamics(::DarwinFactory) = (
    DIC = DIC_geider_light,
    DIN = DIN_geider_light,
    PO4 = PO4_geider_light,
    DOC = DOC_default,
    POC = POC_default,
    DON = DON_default,
    PON = PON_default,
    DOP = DOP_default,
    POP = POP_default,
)

"""Default biogeochemical specification for DARWIN."""
function default_biogeochem_args(::DarwinFactory, ::Type{FT}) where {FT<:AbstractFloat}
    return BiogeochemistrySpecification(
        POC_remineralization = FT(0.1213 / day),
        DOC_remineralization = FT(0.1213 / day),
        PON_remineralization = FT(0.1213 / day),
        DON_remineralization = FT(0.1213 / day),
        POP_remineralization = FT(0.1213 / day),
        DOP_remineralization = FT(0.1213 / day),
        DOM_POM_fractionation = FT(0.45),
        nitrogen_to_carbon = FT(0.15),
        phosphorus_to_carbon = FT(0.009),
    )
end

"""Default interaction builder for DARWIN."""
function default_interactions(::DarwinFactory)
    return InteractionDynamics(_darwin_default_interactions)
end

function _darwin_default_interactions(ctx, args=NamedTuple())
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
