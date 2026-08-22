"""Parameter definitions for the NiPiZD model.

This file defines a single source of truth for:
- parameter definitions and semantic provisions
- constructor-time default values (via `DefaultProvider` entries)

Numeric defaults are evaluated on the host during construction and later moved to the
target architecture with `Adapt`.

Interaction matrices (`palatability_matrix`, `assimilation_matrix`) use derived
defaults whose dependencies are declared beside their parameter definitions.
"""

import ...Parameters:
    parameter_definitions,
    ParameterDefinition,
    ParameterProvision,
    ConstantDefault,
    DerivedDefault,
    DiameterIndexedVectorDefault

using ...Library.Allometry: AllometricParam, PowerLaw

using ...Configuration: PalatabilityAllometric, AssimilationBinary

function parameter_definitions(::NiPiZDFamily)
    detritus_remin = 0.1213 / 86400

    return (
        ParameterDefinition(
            :detritus_remineralization,
            ConstantDefault(detritus_remin);
            provides=ParameterProvision(:remineralization_D, :rate),
        ),
        ParameterDefinition(
            :linear_detrital_mortality,
            DiameterIndexedVectorDefault(0.0101 / 86400; default=0);
            axes=:plankton,
            provides=ParameterProvision(:linear_mortality_P_to_D, :rate),
        ),
        ParameterDefinition(
            :linear_mortality,
            DiameterIndexedVectorDefault(8e-7; default=0);
            axes=:plankton,
            provides=(
                ParameterProvision(:linear_mortality_P_to_N, :rate),
                ParameterProvision(:linear_mortality_Z_to_N, :rate),
            ),
        ),
        ParameterDefinition(
            :quadratic_mortality,
            DiameterIndexedVectorDefault(1e-6; default=0);
            axes=:plankton,
            provides=ParameterProvision(:quadratic_mortality_Z_to_D, :rate),
        ),
        ParameterDefinition(
            :maximum_growth_rate,
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15);
                default=0,
            );
            axes=:plankton,
            provides=ParameterProvision(:growth_P, :maximum_rate),
        ),
        ParameterDefinition(
            :nutrient_half_saturation,
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27);
                default=0,
            );
            axes=:plankton,
            provides=ParameterProvision(:growth_P, :K),
        ),
        ParameterDefinition(
            :alpha,
            DiameterIndexedVectorDefault(
                0.1953 / 86400; default=0
            );
            axes=:plankton,
            provides=ParameterProvision(:growth_P, :alpha),
        ),
        ParameterDefinition(
            :maximum_predation_rate,
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16);
                default=0,
            );
            axes=:plankton,
            provides=ParameterProvision(:grazing_Z_on_P, :maximum_rate),
        ),
        ParameterDefinition(
            :holling_half_saturation,
            DiameterIndexedVectorDefault(5.0; default=0);
            axes=:plankton,
            provides=ParameterProvision(:grazing_Z_on_P, :half_saturation),
        ),
        ParameterDefinition(
            :palatability_matrix,
            DerivedDefault(
                PalatabilityAllometric();
                deps=(
                    :optimum_predator_prey_ratio,
                    :specificity,
                    :protection,
                ),
            );
            axes=(:consumer, :prey),
            provides=ParameterProvision(:grazing_Z_on_P, :palatability),
        ),
        ParameterDefinition(
            :assimilation_matrix,
            DerivedDefault(
                AssimilationBinary(); deps=(:assimilation_efficiency,)
            );
            axes=(:consumer, :prey),
            provides=ParameterProvision(:grazing_Z_on_P, :assimilation),
        ),
        ParameterDefinition(
            :optimum_predator_prey_ratio,
            DiameterIndexedVectorDefault(10.0; default=0);
            axes=:plankton,
        ),
        ParameterDefinition(
            :specificity,
            DiameterIndexedVectorDefault(0.3; default=0);
            axes=:plankton,
        ),
        ParameterDefinition(
            :protection,
            DiameterIndexedVectorDefault(0.0; default=1.0);
            axes=:plankton,
        ),
        ParameterDefinition(
            :assimilation_efficiency,
            DiameterIndexedVectorDefault(0.32; default=0);
            axes=:plankton,
        ),
    )
end
