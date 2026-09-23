import ...Parameters:
    parameter_definitions,
    Parameter,
    ConstructionParameter,
    DerivedDefault,
    DiameterIndexedVectorDefault,
    ConsumerResourceFromConsumer

using ...Library.Allometry: AllometricParam, PowerLaw
using ...Parameters: AllometricPalatability

"""LOBSTER3-like defaults expressed through Agate size-trait machinery."""
function parameter_definitions(::FrankenLOBSTERFamily)
    day = 86400

    # For the canonical P diameters (< 3 um), the supplied LOBSTER3 implementation uses
    # mu = 1.2066 * V^0.28 / day. Its nitrate half-saturation construction reduces to
    # 0.028154 * V^0.65. The supplied Smith/analytical-light slope is 0.1953 / day.
    maximum_growth = AllometricParam(
        PowerLaw(); prefactor=1.2066 / day, exponent=0.28
    )
    nitrate_half_saturation = AllometricParam(
        PowerLaw(); prefactor=0.028154, exponent=0.65
    )
    ammonium_half_saturation = AllometricParam(
        PowerLaw(); prefactor=0.5 * 0.028154, exponent=0.65
    )

    # Supplied LOBSTER3 heterotroph coefficients (Follett/Zakem/DARWIN family):
    # mu_max = 1.836 * V^0.28 / day and
    # K_DOM = k * mu_max * Qmin / Vmax = 0.04284 * V^0.65.
    bacterial_maximum_uptake = AllometricParam(
        PowerLaw(); prefactor=1.836 / day, exponent=0.28
    )
    bacterial_dom_half_saturation = AllometricParam(
        PowerLaw(); prefactor=0.04284, exponent=0.65
    )

    return (
        maximum_growth_rate=Parameter(
            DiameterIndexedVectorDefault(maximum_growth; default=0)
        ),
        nitrate_half_saturation=Parameter(
            DiameterIndexedVectorDefault(nitrate_half_saturation; default=0)
        ),
        ammonium_half_saturation=Parameter(
            DiameterIndexedVectorDefault(ammonium_half_saturation; default=0)
        ),
        iron_half_saturation=Parameter(2e-4),
        ammonium_inhibition=Parameter(3.0),
        temperature_q10=Parameter(1.88),
        reference_temperature=Parameter(20.0),
        alpha=Parameter(DiameterIndexedVectorDefault(0.1953 / day; default=0)),
        phytoplankton_exudation_fraction=Parameter(0.05),
        ammonium_fraction_of_exudate=Parameter(0.75),
        phytoplankton_mortality_rate=Parameter(5.8e-7),
        zooplankton_excretion_rate=Parameter(5.8e-7),
        ammonium_fraction_of_zooplankton_excretion=Parameter(0.5),
        zooplankton_mortality_rate=Parameter(2.31e-6),
        bacterial_maximum_uptake_rate=Parameter(
            DiameterIndexedVectorDefault(bacterial_maximum_uptake; default=0)
        ),
        bacterial_dom_half_saturation=Parameter(
            DerivedDefault(
                ConsumerResourceFromConsumer(); deps=(:bacterial_dom_affinity_trait,)
            )
        ),
        bacterial_substrate_preference=Parameter(1.0),
        bacterial_assimilation=Parameter(0.1),
        bacterioplankton_mortality_rate=Parameter(5.8e-7),
        maximum_predation_rate=Parameter(
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=15.9 / day, exponent=-0.16);
                default=0,
            )
        ),
        grazing_half_saturation=Parameter(1.0),
        palatability_matrix=Parameter(
            DerivedDefault(
                AllometricPalatability();
                deps=(:optimum_predator_prey_ratio, :specificity),
            )
        ),
        assimilation_matrix=Parameter(0.7),
        optimum_predator_prey_ratio=ConstructionParameter(
            DiameterIndexedVectorDefault(10.0; default=0); axes=:plankton
        ),
        specificity=ConstructionParameter(
            DiameterIndexedVectorDefault(0.3; default=0); axes=:plankton
        ),
        bacterial_dom_affinity_trait=ConstructionParameter(
            DiameterIndexedVectorDefault(bacterial_dom_half_saturation; default=0);
            axes=:plankton,
        ),
    )
end
