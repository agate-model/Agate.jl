using ...ModelFamilies: AbstractModelFamily
using ...Components: Plankton, Pool
using ...Processes:
    Growth,
    Light,
    NutrientResponse,
    Consumption,
    Mortality,
    Products,
    Smith,
    Monod,
    PreferentialGrazing,
    HeterotrophicConsumption,
    QuadraticMortality

import ...ModelFamilies: default_components, default_processes, definition_version
import ...Construction: family_id, registered_family

"""Registered family for the Agate-generated FrankenLOBSTER living community."""
struct FrankenLOBSTERFamily <: AbstractModelFamily end

family_id(::FrankenLOBSTERFamily) = :FrankenLOBSTER
registered_family(::Val{:FrankenLOBSTER}) = FrankenLOBSTERFamily()
definition_version(::FrankenLOBSTERFamily)::VersionNumber = v"0.7.0"

"""LOBSTER3-like default living-community size structure."""
const DEFAULT_SIZE_STRUCTURE = (
    phytoplankton=(P=(n=2, min_esd=0.6, max_esd=1.2, spacing=:linear),),
    zooplankton=(Z=(n=2, min_esd=6.0, max_esd=12.0, spacing=:linear),),
    bacterioplankton=(H=(n=1, min_esd=0.6, max_esd=0.6, spacing=:linear),),
)

# NO3, NH4, and DOM are OceanBioME-owned state used by the compiled living-community
# equations. Waste pools are exchange accumulators reported through NPD hooks rather than
# prognostic fields owned by Agate.
const FRANKENLOBSTER_COMPONENTS = (
    NO₃=Pool(:nitrogen),
    NH₄=Pool(:nitrogen),
    DOM=Pool(:nitrogen),
    solid_waste=Pool(:nitrogen),
    inorganic_waste=Pool(:nitrogen),
    P=Plankton(;
        states=(nitrogen=:nitrogen,),
        reference_state=:nitrogen,
        size_structure=DEFAULT_SIZE_STRUCTURE.phytoplankton.P,
    ),
    Z=Plankton(;
        states=(nitrogen=:nitrogen,),
        reference_state=:nitrogen,
        size_structure=DEFAULT_SIZE_STRUCTURE.zooplankton.Z,
    ),
    H=Plankton(;
        states=(nitrogen=:nitrogen,),
        reference_state=:nitrogen,
        size_structure=DEFAULT_SIZE_STRUCTURE.bacterioplankton.H,
    ),
)

default_components(::FrankenLOBSTERFamily) = FRANKENLOBSTER_COMPONENTS

const FRANKENLOBSTER_PROCESSES = (
    # v0.15 deliberately uses nitrate-only phytoplankton growth. Shared-capacity NO3/NH4
    # acquisition is deferred to a follow-up rather than giving the two N sources separate
    # maximum growth capacities.
    nitrate_growth_P=Growth(;
        plankton=:P,
        reference_resource=:NO₃,
        bindings=(maximum_rate=:maximum_growth_rate,),
        factors=(
            light=Light(Smith(); driver=:PAR),
            nutrient=NutrientResponse(
                Monod();
                resource=:NO₃,
                bindings=(half_saturation=:nitrate_half_saturation,),
            ),
        ),
    ),
    consumption_H_on_DOM=Consumption(
        HeterotrophicConsumption();
        consumers=:H,
        resources=:DOM,
        bindings=(
            maximum_rate=:bacterial_maximum_uptake_rate,
            half_saturation=:bacterial_dom_half_saturation,
            substrate_preference=:bacterial_substrate_preference,
            assimilation=:bacterial_assimilation,
        ),
        unassimilated_products=:inorganic_waste,
    ),
    grazing_Z_on_living=Consumption(
        PreferentialGrazing();
        consumers=:Z,
        resources=(:P, :H),
        bindings=(
            maximum_rate=:maximum_predation_rate,
            half_saturation=:grazing_half_saturation,
            palatability=:palatability_matrix,
            assimilation=:assimilation_matrix,
        ),
        unassimilated_products=:solid_waste,
    ),
    mortality_P=Mortality(
        QuadraticMortality();
        plankton=:P,
        bindings=(rate=:phytoplankton_mortality_rate,),
        products=Products(:solid_waste),
    ),
    mortality_Z=Mortality(
        QuadraticMortality();
        plankton=:Z,
        bindings=(rate=:zooplankton_mortality_rate,),
        products=Products(:solid_waste),
    ),
    mortality_H=Mortality(
        QuadraticMortality();
        plankton=:H,
        bindings=(rate=:bacterioplankton_mortality_rate,),
        products=Products(:solid_waste),
    ),
)

default_processes(::FrankenLOBSTERFamily) = FRANKENLOBSTER_PROCESSES
