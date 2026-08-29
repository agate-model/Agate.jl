"""Parameter definitions and construction-time defaults for NiPiZD.

Interaction matrices use `DerivedDefault` providers whose dependencies are declared beside
their parameter definitions.
"""

import ...Parameters:
    parameter_definitions,
    Parameter,
    ConstructionParameter,
    DerivedDefault,
    DiameterIndexedVectorDefault

using ...Library.Allometry: AllometricParam, PowerLaw

using ...Configuration: AllometricPalatability, ConsumerAssimilation

function parameter_definitions(::NiPiZDFamily)
    detritus_remin = 0.1213 / 86400

    return (
        detritus_remineralization=Parameter(detritus_remin),
        mortality_export_fraction=Parameter(0.2),
        linear_mortality=Parameter(
            DiameterIndexedVectorDefault(8e-7; default=0)
        ),
        quadratic_mortality=Parameter(
            DiameterIndexedVectorDefault(1e-6; default=0)
        ),
        maximum_growth_rate=Parameter(
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15);
                default=0,
            )
        ),
        nutrient_half_saturation=Parameter(
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27);
                default=0,
            )
        ),
        alpha=Parameter(
            DiameterIndexedVectorDefault(0.1953 / 86400; default=0)
        ),
        maximum_predation_rate=Parameter(
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16);
                default=0,
            )
        ),
        holling_half_saturation=Parameter(
            DiameterIndexedVectorDefault(5.0; default=0)
        ),
        palatability_matrix=Parameter(
            DerivedDefault(
                AllometricPalatability();
                deps=(
                    :optimum_predator_prey_ratio,
                    :specificity,
                    :protection,
                ),
            )
        ),
        assimilation_matrix=Parameter(
            DerivedDefault(
                ConsumerAssimilation(); deps=(:assimilation_efficiency,)
            )
        ),
        optimum_predator_prey_ratio=ConstructionParameter(
            DiameterIndexedVectorDefault(10.0; default=0);
            axes=:plankton,
        ),
        specificity=ConstructionParameter(
            DiameterIndexedVectorDefault(0.3; default=0);
            axes=:plankton,
        ),
        protection=ConstructionParameter(
            DiameterIndexedVectorDefault(0.0; default=1.0);
            axes=:plankton,
        ),
        assimilation_efficiency=ConstructionParameter(
            DiameterIndexedVectorDefault(0.32; default=0);
            axes=:plankton,
        ),
    )
end
