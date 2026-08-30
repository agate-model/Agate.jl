using ...ModelFamilies: AbstractModelFamily
using ...Components: Plankton, Pool
using ...Processes:
    Growth, Light, NutrientResponse, Consumption, Mortality, Products, Remineralization,
    Smith, Monod, PreferentialGrazing, LinearMortality, QuadraticMortality,
    LinearRemineralization

import ...ModelFamilies: default_components, default_processes, definition_version
import ...Construction: family_id, registered_family

"""Registered family for the size-structured NiPiZD model."""
struct NiPiZDFamily <: AbstractModelFamily end

family_id(::NiPiZDFamily) = :NiPiZD
registered_family(::Val{:NiPiZD}) = NiPiZDFamily()
definition_version(::NiPiZDFamily)::VersionNumber = v"0.1.0"

const DEFAULT_SIZE_STRUCTURE = (
    phytoplankton=(P=(n=2, min_esd=2, max_esd=10, spacing=:log),),
    zooplankton=(Z=(n=2, min_esd=20, max_esd=100, spacing=:linear),),
)

const NIPIZD_COMPONENTS = (
    P=Plankton(;
        states=:nitrogen,
        reference_state=:nitrogen,
        size_structure=DEFAULT_SIZE_STRUCTURE.phytoplankton.P,
    ),
    Z=Plankton(;
        states=:nitrogen,
        reference_state=:nitrogen,
        size_structure=DEFAULT_SIZE_STRUCTURE.zooplankton.Z,
    ),
    N=Pool(:nitrogen),
    D=Pool(:nitrogen),
)

"""Canonical logical components for NiPiZD."""
default_components(::NiPiZDFamily) = NIPIZD_COMPONENTS

const NIPIZD_PROCESSES = (
    growth_P=Growth(;
        plankton=:P,
        reference_resource=:N,
        bindings=(maximum_rate=:maximum_growth_rate,),
        factors=(
            light=Light(Smith(); driver=:PAR),
            nutrients=NutrientResponse(
                Monod(); resource=:N, bindings=(half_saturation=:nutrient_half_saturation,)
            ),
        ),
    ),
    grazing_Z_on_P=Consumption(
        PreferentialGrazing();
        consumers=:Z,
        resources=:P,
        bindings=(
            maximum_rate=:maximum_predation_rate,
            half_saturation=:holling_half_saturation,
            palatability=:palatability_matrix,
            assimilation=:assimilation_matrix,
        ),
        unassimilated_products=:D,
    ),
    linear_mortality_P=Mortality(
        LinearMortality();
        plankton=:P,
        bindings=(rate=:linear_mortality,),
        products=Products(
            (detritus=:D, nutrient=:N);
            fractions=(nutrient=:mortality_export_fraction,),
        ),
    ),
    linear_mortality_Z=Mortality(
        LinearMortality();
        plankton=:Z,
        bindings=(rate=:linear_mortality,),
        products=Products(
            (detritus=:D, nutrient=:N);
            fractions=(nutrient=:mortality_export_fraction,),
        ),
    ),
    quadratic_mortality_Z=Mortality(
        QuadraticMortality();
        plankton=:Z,
        bindings=(rate=:quadratic_mortality,),
        products=Products(
            (detritus=:D, nutrient=:N);
            fractions=(nutrient=:mortality_export_fraction,),
        ),
    ),
    remineralization_D=Remineralization(
        LinearRemineralization();
        sources=:D,
        destination=:N,
        bindings=(rate=:detritus_remineralization,),
    ),
)

"""Canonical named scientific processes for NiPiZD."""
default_processes(::NiPiZDFamily) = NIPIZD_PROCESSES
