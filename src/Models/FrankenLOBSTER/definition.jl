using ...ModelFamilies: AbstractModelFamily
using ...Components: Plankton, Pool
using ...Processes:
    Growth,
    NutrientResponse,
    Consumption,
    Mortality,
    Products,
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
definition_version(::FrankenLOBSTERFamily)::VersionNumber = v"0.4.0"

"""LOBSTER3-like default living-community size structure."""
const DEFAULT_SIZE_STRUCTURE = (
    phytoplankton=(P=(n=2, min_esd=0.6, max_esd=1.2, spacing=:linear),),
    zooplankton=(Z=(n=2, min_esd=6.0, max_esd=12.0, spacing=:linear),),
    bacterioplankton=(B=(n=1, min_esd=0.6, max_esd=0.6, spacing=:linear),),
)

# NO3, NH4, and DOM are OceanBioME-owned state in FrankenLOBSTER. They are represented
# here so Agate processes can use the same named resource identities when compiling the
# living-community equations. `solid_waste` and `inorganic_waste` are exchange
# accumulators rather than fields; their compiled tendencies are reported to NPD.
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
    B=Plankton(;
        states=(nitrogen=:nitrogen,),
        reference_state=:nitrogen,
        size_structure=DEFAULT_SIZE_STRUCTURE.bacterioplankton.B,
    ),
)

"""Canonical logical components for FrankenLOBSTER."""
default_components(::FrankenLOBSTERFamily) = FRANKENLOBSTER_COMPONENTS

const _LIGHT_FACTOR = SaturatingLight()

# Splitting nitrate and ammonium growth into two ordinary Growth processes keeps material
# transfer explicit: each nutrient is removed by exactly the flux that enters phytoplankton.
# Nitrate alone carries the standard LOBSTER ammonium-inhibition factor.
const FRANKENLOBSTER_PROCESSES = (
    nitrate_growth_P=Growth(;
        plankton=:P,
        reference_resource=:NO₃,
        bindings=(maximum_rate=:maximum_growth_rate,),
        factors=(
            light=_LIGHT_FACTOR,
            nutrient=NutrientResponse(
                Monod();
                resource=:NO₃,
                bindings=(half_saturation=:nitrate_half_saturation,),
            ),
            ammonium_inhibition=AmmoniumInhibition(),
        ),
    ),
    ammonium_growth_P=Growth(;
        plankton=:P,
        reference_resource=:NH₄,
        bindings=(maximum_rate=:maximum_growth_rate,),
        factors=(
            light=_LIGHT_FACTOR,
            nutrient=NutrientResponse(
                Monod();
                resource=:NH₄,
                bindings=(half_saturation=:ammonium_half_saturation,),
            ),
        ),
    ),
    consumption_B_on_DOM=Consumption(
        HeterotrophicConsumption();
        consumers=:B,
        resources=:DOM,
        bindings=(
            maximum_rate=:bacterial_maximum_uptake_rate,
            half_saturation=:bacterial_dom_half_saturation,
            substrate_preference=:bacterial_substrate_preference,
            assimilation=:bacterial_assimilation,
        ),
        unassimilated_products=:inorganic_waste,
    ),
    # One grazing process shares each zooplankton ingestion capacity across all living prey.
    grazing_Z_on_living=Consumption(
        PreferentialGrazing();
        consumers=:Z,
        resources=(:P, :B),
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
    mortality_B=Mortality(
        QuadraticMortality();
        plankton=:B,
        bindings=(rate=:bacterioplankton_mortality_rate,),
        products=Products(:solid_waste),
    ),
)

"""Canonical named scientific processes for FrankenLOBSTER."""
default_processes(::FrankenLOBSTERFamily) = FRANKENLOBSTER_PROCESSES
