import ...Parameters:
    parameter_definitions, Parameter, ConstructionParameter, DerivedDefault,
    DiameterIndexedVectorDefault, ConsumerResourceFromConsumer

using ...Library.Allometry: AllometricParam, PowerLaw
using ...Parameters: AllometricPalatability

const FRANKENLOBSTER_CHLOROPHYLL_RATIO = 1.31
const FRANKENLOBSTER_CARBON_RATIO = 6.56
const FRANKENLOBSTER_CALCIUM_CARBONATE_RAIN_RATIO = 0.1
const FRANKENLOBSTER_ZOOPLANKTON_CALCIUM_CARBONATE_DISSOLUTION = 0.3

"""LOBSTER3-like defaults expressed through Agate size-trait machinery."""
function parameter_definitions(::FrankenLOBSTERFamily)
    day = 86400
    law(prefactor, exponent) = AllometricParam(PowerLaw(); prefactor, exponent)
    diameter_default(value) = DiameterIndexedVectorDefault(value; default=0)

    # FrankenLOBSTER keeps size-dependent traits while using LOBSTER light and N responses.
    maximum_growth = law(1.2066 / day, 0.28)
    nitrate_affinity = law(0.028154, 0.65)
    ammonia_affinity = law(0.5 * 0.028154, 0.65)
    bacterial_uptake = law(1.836 / day, 0.28)
    bacterial_affinity = law(0.04284, 0.65)

    return (
        maximum_growth_rate=Parameter(diameter_default(maximum_growth)),
        nitrate_half_saturation=Parameter(diameter_default(nitrate_affinity)),
        ammonia_half_saturation=Parameter(diameter_default(ammonia_affinity)),
        nitrate_ammonia_inhibition=Parameter(3.0),
        light_half_saturation=Parameter(33.0),
        temperature_q10=Parameter(1.88),
        reference_temperature=Parameter(20.0),
        phytoplankton_exudation_fraction=Parameter(0.05),
        ammonium_fraction_of_exudate=Parameter(0.75),
        phytoplankton_mortality_rate=Parameter(5.8e-7),
        zooplankton_excretion_rate=Parameter(5.8e-7),
        ammonium_fraction_of_zooplankton_excretion=Parameter(0.5),
        zooplankton_mortality_rate=Parameter(2.31e-6),
        bacterial_maximum_uptake_rate=Parameter(diameter_default(bacterial_uptake)),
        bacterial_dom_half_saturation=Parameter(DerivedDefault(
            ConsumerResourceFromConsumer(); deps=(:bacterial_dom_affinity_trait,)
        )),
        bacterial_substrate_preference=Parameter(1.0),
        bacterial_assimilation=Parameter(0.1),
        bacterioplankton_mortality_rate=Parameter(5.8e-7),
        maximum_predation_rate=Parameter(diameter_default(law(15.9 / day, -0.16))),
        grazing_half_saturation=Parameter(1.0),
        palatability_matrix=Parameter(DerivedDefault(
            AllometricPalatability(); deps=(:optimum_predator_prey_ratio, :specificity)
        )),
        assimilation_matrix=Parameter(0.7),
        optimum_predator_prey_ratio=ConstructionParameter(
            diameter_default(10.0); axes=:plankton
        ),
        specificity=ConstructionParameter(diameter_default(0.3); axes=:plankton),
        bacterial_dom_affinity_trait=ConstructionParameter(
            diameter_default(bacterial_affinity); axes=:plankton
        ),
    )
end
