using ...ModelFamilies: AbstractModelFamily
using ...Configuration: Population, Pool
using ...Processes:
    Growth, Light, NutrientResponse, Consumption, Mortality, ProductRouting, Remineralization

import ...ModelFamilies: default_components, default_processes
import ...Construction: family_id, registered_family

"""Registered family for the size-structured NiPiZD model."""
struct NiPiZDFamily <: AbstractModelFamily end

family_id(::NiPiZDFamily) = :NiPiZD
registered_family(::Val{:NiPiZD}) = NiPiZDFamily()

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
default_components(::NiPiZDFamily) = NIPIZD_COMPONENTS

const NIPIZD_PROCESSES = (
    growth_P=Growth(;
        population=:P,
        factors=(
            light=Light(:smith; driver=:PAR),
            nutrients=NutrientResponse(:monod; resource=:N),
        ),
    ),
    grazing_Z_on_P=Consumption(
        :idealized; consumer=:Z, resource=:P, unassimilated_destination=:D
    ),
    linear_mortality_P_to_N=Mortality(
        :linear;
        population=:P,
        routing=ProductRouting(:direct; destination=:N),
    ),
    linear_mortality_P_to_D=Mortality(
        :linear;
        population=:P,
        routing=ProductRouting(:direct; destination=:D),
    ),
    linear_mortality_Z_to_N=Mortality(
        :linear;
        population=:Z,
        routing=ProductRouting(:direct; destination=:N),
    ),
    quadratic_mortality_Z_to_D=Mortality(
        :quadratic;
        population=:Z,
        routing=ProductRouting(:direct; destination=:D),
    ),
    remineralization_D=Remineralization(:linear; source=:D, destination=:N),
)

"""Canonical named scientific processes for NiPiZD."""
default_processes(::NiPiZDFamily) = NIPIZD_PROCESSES
