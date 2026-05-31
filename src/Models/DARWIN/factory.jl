using ...Factories: AbstractBGCFactory
using ...Configuration: PFTSpecification

# NOTE: Numeric parameter defaults are declared in `Models/DARWIN/parameters.jl`.

import ...Factories:
    default_plankton_dynamics, default_community, default_biogeochem_dynamics

using ...Tendencies:
    TendencyConfig,
    nutrient_coupling,
    phytoplankton_tendency,
    zooplankton_tendency,
    inorganic_tendency,
    organic_matter_tendency

"""Factory for the simplified DARWIN-like elemental cycling model."""
struct DarwinFactory <: AbstractBGCFactory end

const DARWIN_TENDENCIES = TendencyConfig(;
    growth=:geider,
    organic_cycling=:dom_pom,
    zooplankton=:preferential_grazing,
    nutrient_limitation=:liebig,
    nutrients=(
        nutrient_coupling(
            :DIN,
            :half_saturation_DIN;
            stoichiometry=:nitrogen_to_carbon,
            remineralization=((:DON, :DON_remineralization), (:PON, :PON_remineralization)),
        ),
        nutrient_coupling(
            :PO4,
            :half_saturation_PO4;
            stoichiometry=:phosphorus_to_carbon,
            remineralization=((:DOP, :DOP_remineralization), (:POP, :POP_remineralization)),
        ),
    ),
)

phytoplankton_darwin(plankton_idx::Int) = phytoplankton_tendency(
    DARWIN_TENDENCIES; plankton_idx
)

zooplankton_darwin(plankton_idx::Int) = zooplankton_tendency(
    DARWIN_TENDENCIES; plankton_idx
)

"""Default plankton dynamics for DARWIN."""
function default_plankton_dynamics(::DarwinFactory)
    (Z=zooplankton_darwin, P=phytoplankton_darwin)
end

"""Default structural parameter arguments for DARWIN.

Returns a `NamedTuple` mapping group prefix => group specification.

Ordering is significant; the default keeps the historical `Z`-then-`P` ordering.
"""
function default_community(::DarwinFactory)
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

"""Default non-plankton tracer dynamics for DARWIN."""
function default_biogeochem_dynamics(::DarwinFactory)
    (
        DIC=() -> inorganic_tendency(
            DARWIN_TENDENCIES;
            target=:DIC,
            remineralization=((:DOC, :DOC_remineralization), (:POC, :POC_remineralization)),
            stoichiometry=:one,
        ),
        DIN=() -> inorganic_tendency(DARWIN_TENDENCIES; target=:DIN),
        PO4=() -> inorganic_tendency(DARWIN_TENDENCIES; target=:PO4),
        DOC=() -> organic_matter_tendency(
            DARWIN_TENDENCIES; target=:DOC, remineralization=:DOC_remineralization, fraction=:DOM
        ),
        POC=() -> organic_matter_tendency(
            DARWIN_TENDENCIES; target=:POC, remineralization=:POC_remineralization, fraction=:POM
        ),
        DON=() -> organic_matter_tendency(
            DARWIN_TENDENCIES;
            target=:DON,
            remineralization=:DON_remineralization,
            fraction=:DOM,
            stoichiometry=:nitrogen_to_carbon,
        ),
        PON=() -> organic_matter_tendency(
            DARWIN_TENDENCIES;
            target=:PON,
            remineralization=:PON_remineralization,
            fraction=:POM,
            stoichiometry=:nitrogen_to_carbon,
        ),
        DOP=() -> organic_matter_tendency(
            DARWIN_TENDENCIES;
            target=:DOP,
            remineralization=:DOP_remineralization,
            fraction=:DOM,
            stoichiometry=:phosphorus_to_carbon,
        ),
        POP=() -> organic_matter_tendency(
            DARWIN_TENDENCIES;
            target=:POP,
            remineralization=:POP_remineralization,
            fraction=:POM,
            stoichiometry=:phosphorus_to_carbon,
        ),
    )
end
