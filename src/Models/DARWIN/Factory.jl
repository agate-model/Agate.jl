"""DARWIN defaults expressed as a model-agnostic factory.

This file defines `DarwinFactory` and the default inputs used by
`Agate.Constructor.construct_factory(factory; ...)`.

**Note:** All parameter defaults live exclusively in the DARWIN default parameter
generator (see `Models/DARWIN/Parameters.jl`). This factory provides only
structural defaults (community sizes/diameters) and default dynamics functions.
"""

using ...Utils: AbstractBGCFactory
using ...Utils.Specifications: PFTSpecification
using ...Utils: DiameterRangeSpecification

# NOTE: Defaults are provided by `Constructor.default_parameters` (see `Models/DARWIN/Parameters.jl`).

import ...FactoryInterface:
    default_plankton_dynamics,
    default_community,
    default_biogeochem_dynamics,
    factory_groups,
    consumer_groups,
    prey_groups

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

"""Return the fixed group set for DARWIN.

The ordering matches `default_community` (`Z` then `P`).
"""
factory_groups(::DarwinFactory) = (:Z, :P)

"""Plankton groups that consume in predator-by-prey matrices."""
consumer_groups(::DarwinFactory) = (:Z,)

"""Plankton groups that may be eaten in predator-by-prey matrices."""
prey_groups(::DarwinFactory) = (:P,)

"""Default plankton dynamics for DARWIN."""
default_plankton_dynamics(::DarwinFactory) =
    (Z=zooplankton_default, P=phytoplankton_growth_two_nutrients_geider_light)

"""Default structural parameter arguments for DARWIN.

Returns a `NamedTuple` mapping group prefix => group specification.

Ordering is significant; the default keeps the historical `Z`-then-`P` ordering.
"""
function default_community(::DarwinFactory)
    # Structural defaults only (sizes/diameters). No parameter defaults.
    empty_pft = PFTSpecification()
    return (
        Z=(;
            n=2,
            diameters=DiameterRangeSpecification(20, 100, :linear_splitting),
            pft=empty_pft,
        ),
        P=(;
            n=2, diameters=DiameterRangeSpecification(2, 10, :log_splitting), pft=empty_pft
        ),
    )
end

"""Default non-plankton tracer dynamics for DARWIN."""
default_biogeochem_dynamics(::DarwinFactory) = (
    DIC=DIC_geider_light,
    DIN=DIN_geider_light,
    PO4=PO4_geider_light,
    DOC=DOC_default,
    POC=POC_default,
    DON=DON_default,
    PON=PON_default,
    DOP=DOP_default,
    POP=POP_default,
)
