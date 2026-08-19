using ...Factories: AbstractBGCFactory
using ...Configuration: PFTSpecification, Population, Pool
using ...Processes:
    Growth, NutrientResponse, Grazing, Mortality, ProductRouting, Remineralization

# NOTE: Numeric parameter defaults are declared alongside parameter metadata in
# `parameter_definitions(::NiPiZDFactory)` (see `Models/NiPiZD/parameters.jl`).

import ...Factories:
    default_components,
    default_processes,
    default_plankton_dynamics,
    default_community,
    default_biogeochem_dynamics

import ...Construction: recipe_family, recipe_factory

using ...Tendencies:
    TendencyConfig,
    nutrient_coupling,
    phytoplankton_tendency,
    zooplankton_tendency,
    inorganic_tendency,
    detritus_tendency

"""Factory for the size-structured NiPiZD model."""
struct NiPiZDFactory <: AbstractBGCFactory end

recipe_family(::NiPiZDFactory) = :NiPiZD
recipe_factory(::Val{:NiPiZD}) = NiPiZDFactory()

const DEFAULT_SIZE_STRUCTURE = (
    phytoplankton=(P=(n=2, min_esd=2, max_esd=10, splitting=:log_splitting),),
    zooplankton=(Z=(n=2, min_esd=20, max_esd=100, splitting=:linear_splitting),),
)

const NIPIZD_COMPONENTS = (
    P=Population(;
        currency=:nitrogen, size_structure=DEFAULT_SIZE_STRUCTURE.phytoplankton.P
    ),
    Z=Population(;
        currency=:nitrogen, size_structure=DEFAULT_SIZE_STRUCTURE.zooplankton.Z
    ),
    N=Pool(:nitrogen),
    D=Pool(:nitrogen),
)

"""Canonical logical components for NiPiZD."""
default_components(::NiPiZDFactory) = NIPIZD_COMPONENTS

const NIPIZD_PROCESSES = (
    growth_P=Growth(
        :smith;
        population=:P,
        light=:PAR,
        limitation=NutrientResponse(:monod; resource=:N),
    ),
    grazing_Z_on_P=Grazing(
        :preferential; consumer=:Z, resource=:P, unassimilated_destination=:D
    ),
    linear_mortality_P=Mortality(
        :linear;
        population=:P,
        routing=ProductRouting(:partition; retained=:D, exported=:N),
    ),
    linear_mortality_Z=Mortality(
        :linear;
        population=:Z,
        routing=ProductRouting(:partition; retained=:D, exported=:N),
    ),
    quadratic_mortality_Z=Mortality(
        :quadratic;
        population=:Z,
        routing=ProductRouting(:partition; retained=:D, exported=:N),
    ),
    remineralization_D=Remineralization(:linear; source=:D, destination=:N),
)

"""Canonical named scientific processes for NiPiZD."""
default_processes(::NiPiZDFactory) = NIPIZD_PROCESSES

const NIPIZD_TENDENCIES = TendencyConfig(;
    growth=:smith,
    organic_cycling=:simple_detritus,
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

function default_plankton_dynamics(
    ::NiPiZDFactory, community::NamedTuple, ecological_roles::NamedTuple
)
    groups = keys(community)
    phytoplankton = ecological_roles.phytoplankton
    zooplankton = ecological_roles.zooplankton

    values = ntuple(length(groups)) do i
        group = groups[i]
        if group in phytoplankton
            phytoplankton_nipizd
        elseif group in zooplankton
            zooplankton_nipizd
        else
            throw(ArgumentError("recipe group :$group has no NiPiZD ecological role"))
        end
    end
    return NamedTuple{groups}(values)
end

"""Default plankton arguments for NiPiZD.

Returns a `NamedTuple` mapping group prefix => group specification.

Ordering is significant; the default uses `Z`-then-`P` ordering.
"""
function default_community(::NiPiZDFactory)
    # Structural defaults only (sizes/diameters). No parameter defaults.
    empty_pft = PFTSpecification()
    return (
        Z=(; diameters=DEFAULT_SIZE_STRUCTURE.zooplankton.Z, pft=empty_pft),
        P=(; diameters=DEFAULT_SIZE_STRUCTURE.phytoplankton.P, pft=empty_pft),
    )
end

"""Default non-plankton tracer dynamics for NiPiZD."""
function default_biogeochem_dynamics(::NiPiZDFactory)
    (
        N=() -> inorganic_tendency(NIPIZD_TENDENCIES; target=:N),
        D=() -> detritus_tendency(NIPIZD_TENDENCIES),
    )
end
