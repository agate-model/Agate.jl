"""NiPiZD defaults expressed as a model-agnostic factory.

This file defines `NiPiZDFactory` and the default inputs used by
`Agate.Constructor.construct(factory; ...)`.
"""

using Oceananigans.Units

using Agate.Utils: AbstractBGCFactory
using Agate.Utils.Specifications: PFTSpecification, BiogeochemistrySpecification
using Agate.Utils: DiameterRangeSpecification

using Agate.Library.Allometry: AllometricParam, PowerLaw

import Agate.Models: default_plankton_dynamics, default_plankton_args, default_biogeochem_dynamics, default_biogeochem_args
using .Tracers:
    nutrient_default,
    detritus_default,
    phytoplankton_default,
    zooplankton_default


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
    phyto_pft = PFTSpecification(
        maximum_growth_rate = AllometricParam(PowerLaw(); prefactor=FT(2 / day), exponent=FT(-0.15)),
        nutrient_half_saturation = AllometricParam(PowerLaw(); prefactor=FT(0.17), exponent=FT(0.27)),
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

    zoo_pft = PFTSpecification(
        maximum_predation_rate = AllometricParam(PowerLaw(); prefactor=FT(30.84 / day), exponent=FT(-0.16)),
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


