using ...Factories: AbstractBGCFactory
using ...Configuration: PFTSpecification

# NOTE: Numeric parameter defaults are declared alongside parameter metadata in
# `parameter_definitions(::NiPiZDFactory)` (see `Models/NiPiZD/parameters.jl`).

import ...Factories:
    default_plankton_dynamics, default_community, default_biogeochem_dynamics

using ...Tendencies:
    TendencyConfig,
    nutrient_coupling,
    phytoplankton_tendency,
    zooplankton_tendency,
    inorganic_tendency,
    detritus_tendency

"""Factory for the size-structured NiPiZD model."""
struct NiPiZDFactory <: AbstractBGCFactory end

const NIPIZD_TENDENCIES = TendencyConfig(;
    formulation=:smith_detritus,
    zooplankton=:preferential_grazing,
    nutrient_limitation=:liebig,
    nutrients=(
        nutrient_coupling(
            :N,
            :nutrient_half_saturation;
            remineralization=((:D, :detritus_remineralization),),
        ),
    ),
)

phytoplankton_nipizd(plankton_idx::Int) = phytoplankton_tendency(
    NIPIZD_TENDENCIES; plankton_idx
)

zooplankton_nipizd(plankton_idx::Int) = zooplankton_tendency(
    NIPIZD_TENDENCIES; plankton_idx
)

"""Default plankton dynamics for NiPiZD.

Returns a `NamedTuple` mapping group prefix => tracer dynamics builder.
"""
function default_plankton_dynamics(::NiPiZDFactory)
    (Z=zooplankton_nipizd, P=phytoplankton_nipizd)
end

"""Default plankton arguments for NiPiZD.

Returns a `NamedTuple` mapping group prefix => group specification.

Ordering is significant; the default keeps the historical `Z`-then-`P` ordering.
"""
function default_community(::NiPiZDFactory)
    # Structural defaults only (sizes/diameters). No parameter defaults.
    empty_pft = PFTSpecification()
    return (
        Z=(;
            diameters=(n=2, min_esd=20, max_esd=100, splitting=:linear_splitting),
            pft=empty_pft,
        ),
        P=(;
            diameters=(n=2, min_esd=2, max_esd=10, splitting=:log_splitting), pft=empty_pft
        ),
    )
end

"""Default non-plankton tracer dynamics for NiPiZD."""
function default_biogeochem_dynamics(::NiPiZDFactory)
    (
        N=() -> inorganic_tendency(NIPIZD_TENDENCIES; target=:N),
        D=() -> detritus_tendency(NIPIZD_TENDENCIES),
    )
end
