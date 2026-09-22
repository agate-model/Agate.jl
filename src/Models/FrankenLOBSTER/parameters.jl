import ...Parameters:
    parameter_definitions,
    Parameter,
    ConstructionParameter,
    DerivedDefault,
    DiameterIndexedVectorDefault

using ...Library.Allometry: AllometricParam, PowerLaw
using ...Parameters: AllometricPalatability, ConsumerAssimilation

"""LOBSTER3-like defaults expressed through Agate size-trait machinery."""
function parameter_definitions(::FrankenLOBSTERFamily)
    day = 86400

    # For the canonical P diameters (< 3 um), the supplied LOBSTER3 implementation uses
    # mu = 1.2066 * V^0.28 / day. Its nitrate half-saturation construction combines
    # k * mu * Qmin / Vmax, which reduces to 0.028154 * V^0.65 in this size regime.
    # Following the DARWIN formulation used by Zakem et al. (2018), ammonium affinity
    # is twice nitrate affinity, represented as K_NH4 = 0.5 K_NO3 at every cell size.
    maximum_growth = AllometricParam(
        PowerLaw(); prefactor=1.2066 / day, exponent=0.28
    )
    nitrate_half_saturation = AllometricParam(
        PowerLaw(); prefactor=0.028154, exponent=0.65
    )
    ammonium_half_saturation = AllometricParam(
        PowerLaw(); prefactor=0.5 * 0.028154, exponent=0.65
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
        light_half_saturation=Parameter(55.0),
        nitrate_ammonia_inhibition=Parameter(3.0),
        phytoplankton_mortality_rate=Parameter(5.8e-7),
        zooplankton_mortality_rate=Parameter(2.31e-6),
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
                deps=(:optimum_predator_prey_ratio, :specificity, :protection),
            )
        ),
        assimilation_matrix=Parameter(
            DerivedDefault(
                ConsumerAssimilation(); deps=(:assimilation_efficiency,)
            )
        ),
        optimum_predator_prey_ratio=ConstructionParameter(
            DiameterIndexedVectorDefault(10.0; default=0); axes=:plankton
        ),
        specificity=ConstructionParameter(
            DiameterIndexedVectorDefault(0.3; default=0); axes=:plankton
        ),
        protection=ConstructionParameter(
            DiameterIndexedVectorDefault(0.0; default=1.0); axes=:plankton
        ),
        assimilation_efficiency=ConstructionParameter(
            DiameterIndexedVectorDefault(0.7; default=0); axes=:plankton
        ),
    )
end
