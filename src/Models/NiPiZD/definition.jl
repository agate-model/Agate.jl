using ...ModelFamilies: AbstractModelFamily
using ...Configuration: Population, Pool
using ...Processes:
    Growth, Light, NutrientResponse, Consumption, Mortality, ProductRouting, Remineralization,
    Smith, Monod, IdealizedGrazing, LinearMortality, QuadraticMortality,
    LinearRemineralization, DirectRouting

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
    P=Population(:nitrogen; size_structure=DEFAULT_SIZE_STRUCTURE.phytoplankton.P),
    Z=Population(:nitrogen; size_structure=DEFAULT_SIZE_STRUCTURE.zooplankton.Z),
    N=Pool(:nitrogen),
    D=Pool(:nitrogen),
)

"""Canonical logical components for NiPiZD."""
default_components(::NiPiZDFamily) = NIPIZD_COMPONENTS

const NIPIZD_PROCESSES = (
    growth_P=Growth(;
        populations=:P,
        factors=(
            light=Light(
                Smith(); driver=:PAR, bindings=(maximum_rate=:maximum_growth_rate,)
            ),
            nutrients=NutrientResponse(
                Monod(); resource=:N, bindings=(K=:nutrient_half_saturation,)
            ),
        ),
    ),
    grazing_Z_on_P=Consumption(
        IdealizedGrazing();
        consumers=:Z,
        resources=:P,
        bindings=(
            maximum_rate=:maximum_predation_rate,
            half_saturation=:holling_half_saturation,
            palatability=:palatability_matrix,
            assimilation=:assimilation_matrix,
        ),
        routing=ProductRouting(DirectRouting(); destination=:D),
    ),
    linear_mortality_P_to_N=Mortality(
        LinearMortality();
        populations=:P,
        bindings=(rate=:linear_mortality,),
        routing=ProductRouting(DirectRouting(); destination=:N),
    ),
    linear_mortality_P_to_D=Mortality(
        LinearMortality();
        populations=:P,
        bindings=(rate=:linear_detrital_mortality,),
        routing=ProductRouting(DirectRouting(); destination=:D),
    ),
    linear_mortality_Z_to_N=Mortality(
        LinearMortality();
        populations=:Z,
        bindings=(rate=:linear_mortality,),
        routing=ProductRouting(DirectRouting(); destination=:N),
    ),
    quadratic_mortality_Z_to_D=Mortality(
        QuadraticMortality();
        populations=:Z,
        bindings=(rate=:quadratic_mortality,),
        routing=ProductRouting(DirectRouting(); destination=:D),
    ),
    remineralization_D=Remineralization(
        LinearRemineralization();
        sources=:D,
        destinations=:N,
        bindings=(rate=:detritus_remineralization,),
    ),
)

"""Canonical named scientific processes for NiPiZD."""
default_processes(::NiPiZDFamily) = NIPIZD_PROCESSES
