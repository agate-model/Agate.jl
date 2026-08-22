"""Parameter definitions for the NiPiZD model.

This file defines a single source of truth for:
- keyed parameter definitions and constructor-time defaults

Numeric defaults are evaluated on the host during construction and later moved to the
target architecture with `Adapt`.

Interaction matrices (`palatability_matrix`, `assimilation_matrix`) use derived
defaults whose dependencies are declared beside their parameter definitions.
"""

import ...Parameters:
    parameter_definitions,
    Parameter,
    ConstantDefault,
    DerivedDefault,
    DiameterIndexedVectorDefault

using ...Library.Allometry: AllometricParam, PowerLaw

using ...Configuration: PalatabilityAllometric, AssimilationBinary

function parameter_definitions(::NiPiZDFamily)
    detritus_remin = 0.1213 / 86400

    return (
        detritus_remineralization=Parameter(
            ConstantDefault(detritus_remin),
        ),
        linear_detrital_mortality=Parameter(
            DiameterIndexedVectorDefault(0.0101 / 86400; default=0);
            axes=:plankton,
        ),
        linear_mortality=Parameter(
            DiameterIndexedVectorDefault(8e-7; default=0);
            axes=:plankton,
        ),
        quadratic_mortality=Parameter(
            DiameterIndexedVectorDefault(1e-6; default=0);
            axes=:plankton,
        ),
        maximum_growth_rate=Parameter(
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15);
                default=0,
            );
            axes=:plankton,
        ),
        nutrient_half_saturation=Parameter(
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27);
                default=0,
            );
            axes=:plankton,
        ),
        alpha=Parameter(
            DiameterIndexedVectorDefault(0.1953 / 86400; default=0);
            axes=:plankton,
        ),
        maximum_predation_rate=Parameter(
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16);
                default=0,
            );
            axes=:plankton,
        ),
        holling_half_saturation=Parameter(
            DiameterIndexedVectorDefault(5.0; default=0);
            axes=:plankton,
        ),
        palatability_matrix=Parameter(
            DerivedDefault(
                PalatabilityAllometric();
                deps=(
                    :optimum_predator_prey_ratio,
                    :specificity,
                    :protection,
                ),
            );
            axes=(:consumer, :prey),
        ),
        assimilation_matrix=Parameter(
            DerivedDefault(
                AssimilationBinary(); deps=(:assimilation_efficiency,)
            );
            axes=(:consumer, :prey),
        ),
        optimum_predator_prey_ratio=Parameter(
            DiameterIndexedVectorDefault(10.0; default=0);
            axes=:plankton,
        ),
        specificity=Parameter(
            DiameterIndexedVectorDefault(0.3; default=0);
            axes=:plankton,
        ),
        protection=Parameter(
            DiameterIndexedVectorDefault(0.0; default=1.0);
            axes=:plankton,
        ),
        assimilation_efficiency=Parameter(
            DiameterIndexedVectorDefault(0.32; default=0);
            axes=:plankton,
        ),
    )
end
