using ...ModelFamilies: AbstractModelFamily
using ...Components: Plankton, Pool
using ...Processes:
    Growth, Light, Consumption, Mortality, Products, ExponentialSaturation,
    NutrientResponse, Monod, InhibitedMonod, Temperature, Q10, PreferentialGrazing,
    HeterotrophicConsumption, LinearMortality, QuadraticMortality

import ...ModelFamilies: default_components, default_processes, definition_version, plankton_roles
import ...Construction: family_id, registered_family

struct FrankenLOBSTERFamily <: AbstractModelFamily end
family_id(::FrankenLOBSTERFamily) = :FrankenLOBSTER
registered_family(::Val{:FrankenLOBSTER}) = FrankenLOBSTERFamily()
definition_version(::FrankenLOBSTERFamily)::VersionNumber = v"0.14.0"
plankton_roles(::FrankenLOBSTERFamily) = (
    phytoplankton=:P, zooplankton=:Z, bacterioplankton=:H,
)

const DEFAULT_SIZE_STRUCTURE = (
    phytoplankton=(P=(n=2, min_esd=0.6, max_esd=1.2, spacing=:linear),),
    zooplankton=(Z=(n=2, min_esd=6.0, max_esd=12.0, spacing=:linear),),
    bacterioplankton=(H=(n=1, min_esd=0.6, max_esd=0.6, spacing=:linear),),
)

_nitrogen_plankton(size_structure) = Plankton(;
    states=(nitrogen=:nitrogen,), reference_state=:nitrogen, size_structure
)

# NO3, NH4, T, and DOM are external state; waste pools are NPD exchange accumulators.
const FRANKENLOBSTER_COMPONENTS = (
    NO₃=Pool(:nitrogen), NH₄=Pool(:nitrogen), T=Pool(:temperature), DOM=Pool(:nitrogen),
    solid_waste=Pool(:nitrogen), inorganic_waste=Pool(:nitrogen), dissolved_waste=Pool(:nitrogen),
    P=_nitrogen_plankton(DEFAULT_SIZE_STRUCTURE.phytoplankton.P),
    Z=_nitrogen_plankton(DEFAULT_SIZE_STRUCTURE.zooplankton.Z),
    H=_nitrogen_plankton(DEFAULT_SIZE_STRUCTURE.bacterioplankton.H),
)
default_components(::FrankenLOBSTERFamily) = FRANKENLOBSTER_COMPONENTS

const _P_GROWTH_FACTORS = (
    light=Light(ExponentialSaturation(); driver=:PAR, bindings=(half_saturation=:light_half_saturation,)),
    temperature=Temperature(
        Q10(:plankton); component=:T,
        bindings=(q10=:temperature_q10, reference_temperature=:reference_temperature),
    ),
)
const _EXUDATE_PRODUCTS = Products(
    (dissolved=:dissolved_waste, inorganic=:inorganic_waste);
    fractions=(inorganic=:ammonium_fraction_of_exudate,),
)

function _p_growth(resource, response)
    return Growth(;
        plankton=:P, reference_resource=resource,
        bindings=(maximum_rate=:maximum_growth_rate, product_fraction=:phytoplankton_exudation_fraction),
        factors=merge(_P_GROWTH_FACTORS, (nutrients=response,)), products=_EXUDATE_PRODUCTS,
    )
end

_n_mortality(plankton, rate) = Mortality(
    QuadraticMortality(); plankton, bindings=(rate=rate,), products=Products(:solid_waste)
)

const FRANKENLOBSTER_PROCESSES = (
    nitrate_growth_P=_p_growth(:NO₃, NutrientResponse(
        InhibitedMonod(); resource=:NO₃, inhibitor=:NH₄,
        bindings=(half_saturation=:nitrate_half_saturation, inhibition=:nitrate_ammonia_inhibition),
    )),
    ammonia_growth_P=_p_growth(:NH₄, NutrientResponse(
        Monod(); resource=:NH₄, bindings=(half_saturation=:ammonia_half_saturation,)
    )),
    consumption_H_on_DOM=Consumption(
        HeterotrophicConsumption(); consumers=:H, resources=:DOM,
        bindings=(
            maximum_rate=:bacterial_maximum_uptake_rate,
            half_saturation=:bacterial_dom_half_saturation,
            substrate_preference=:bacterial_substrate_preference,
            assimilation=:bacterial_assimilation,
        ),
        unassimilated_products=:inorganic_waste,
    ),
    grazing_Z_on_living=Consumption(
        PreferentialGrazing(); consumers=:Z, resources=(:P, :H),
        bindings=(
            maximum_rate=:maximum_predation_rate, half_saturation=:grazing_half_saturation,
            palatability=:palatability_matrix, assimilation=:assimilation_matrix,
        ),
        unassimilated_products=:solid_waste,
    ),
    excretion_Z=Mortality(
        LinearMortality(); plankton=:Z, bindings=(rate=:zooplankton_excretion_rate,),
        products=Products(
            (dissolved=:dissolved_waste, inorganic=:inorganic_waste);
            fractions=(inorganic=:ammonium_fraction_of_zooplankton_excretion,),
        ),
    ),
    mortality_P=_n_mortality(:P, :phytoplankton_mortality_rate),
    mortality_Z=_n_mortality(:Z, :zooplankton_mortality_rate),
    mortality_H=_n_mortality(:H, :bacterioplankton_mortality_rate),
)

default_processes(::FrankenLOBSTERFamily) = FRANKENLOBSTER_PROCESSES
