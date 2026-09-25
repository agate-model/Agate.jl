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

using ...Parameters: AllometricPalatability, ConsumerAssimilation

function parameter_definitions(::NiPiZDFamily)
    detritus_remin = 0.1213 / 86400
    law(prefactor, exponent) = AllometricParam(PowerLaw(); prefactor, exponent)
    diameter_default(value; default=0) = DiameterIndexedVectorDefault(value; default)

    return (
        detritus_remineralization=Parameter(detritus_remin),
        mortality_export_fraction=Parameter(0.2),
        linear_mortality=Parameter(diameter_default(8e-7)),
        quadratic_mortality=Parameter(diameter_default(1e-6)),
        maximum_growth_rate=Parameter(diameter_default(law(2 / 86400, -0.15))),
        nutrient_half_saturation=Parameter(diameter_default(law(0.17, 0.27))),
        alpha=Parameter(diameter_default(0.1953 / 86400)),
        maximum_predation_rate=Parameter(diameter_default(law(30.84 / 86400, -0.16))),
        holling_half_saturation=Parameter(diameter_default(5.0)),
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
            diameter_default(10.0); axes=:plankton,
        ),
        specificity=ConstructionParameter(
            diameter_default(0.3); axes=:plankton,
        ),
        protection=ConstructionParameter(
            diameter_default(0.0; default=1.0); axes=:plankton,
        ),
        assimilation_efficiency=ConstructionParameter(
            diameter_default(0.32); axes=:plankton,
        ),
    )
end
