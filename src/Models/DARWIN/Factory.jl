"""DARWIN defaults expressed as a model-agnostic factory.

This file defines `DarwinFactory` and the default inputs used by
`Agate.Constructor.construct(factory; ...)`.
"""

using Oceananigans.Units

using Agate.Utils: AbstractBGCFactory
using Agate.Utils.Specifications: PFTSpecification, BiogeochemistrySpecification
using Agate.Utils: DiameterRangeSpecification

using Agate.Library.Allometry: AllometricParam, PowerLaw

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
        maximum_growth_rate = AllometricParam(PowerLaw(); prefactor=FT(2 / day), exponent=FT(-0.15)),
        half_saturation_DIN = AllometricParam(PowerLaw(); prefactor=FT(0.17), exponent=FT(0.27)),
        half_saturation_PO4 = AllometricParam(PowerLaw(); prefactor=FT(0.17), exponent=FT(0.27)),
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

