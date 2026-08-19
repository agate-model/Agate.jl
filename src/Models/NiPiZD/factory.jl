using ...Factories: AbstractBGCFactory
using ...Configuration: PFTSpecification, Population, Pool
using ...Processes:
    Growth, NutrientResponse, Grazing, Mortality, ProductRouting, Remineralization

import ...Factories: default_components, default_processes, default_community
import ...Construction: recipe_family, recipe_factory

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

"""Default population realization for NiPiZD."""
function default_community(::NiPiZDFactory)
    empty_pft = PFTSpecification()
    return (
        Z=(; diameters=DEFAULT_SIZE_STRUCTURE.zooplankton.Z, pft=empty_pft),
        P=(; diameters=DEFAULT_SIZE_STRUCTURE.phytoplankton.P, pft=empty_pft),
    )
end
