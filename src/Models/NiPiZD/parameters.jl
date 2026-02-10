"""Parameter definitions for the NiPiZD model.

This file defines a single source of truth for:
- parameter metadata (`ParameterSpec`)
- constructor-time default values (via `DefaultProvider` entries)

Numeric defaults are evaluated on the host during construction and later moved to the
target architecture with `Adapt`.

Interaction matrices (`palatability_matrix`, `assimilation_matrix`) are derived from
trait vectors using `matrix_definitions`.
"""

import ...Factories:
    parameter_definitions,
    ParameterDefinition,
    ParameterSpec,
    ConstDefault,
    NoDefault,
    FillDefault,
    DiameterIndexedVectorDefault

import ...Configuration: matrix_definitions

using ...Library.Allometry:
    AllometricParam,
    PowerLaw

using ...Configuration: MatrixDefinition, PalatabilityAllometric, AssimilationBinary

function parameter_definitions(::NiPiZDFactory)
    detritus_remin = 0.1213 / 86400

    return (
        ParameterDefinition(
            ParameterSpec(
                :detritus_remineralization,
                :scalar;
                doc="Detritus remineralization rate.",
            ),
            ConstDefault(detritus_remin),
        ),
        ParameterDefinition(
            ParameterSpec(
                :mortality_export_fraction,
                :scalar;
                doc="Fraction of mortality routed to detritus export.",
            ),
            ConstDefault(0.2),
        ),
        ParameterDefinition(
            ParameterSpec(
                :linear_mortality,
                :vector;
                doc="Linear mortality coefficient per plankton class.",
            ),
            FillDefault(8e-7),
        ),
        ParameterDefinition(
            ParameterSpec(
                :quadratic_mortality,
                :vector;
                doc="Quadratic mortality coefficient per plankton class.",
            ),
            DiameterIndexedVectorDefault(1e-6, :default_consumer_indices; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :maximum_growth_rate,
                :vector;
                doc="Maximum phytoplankton growth rate per plankton class.",
            ),
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15),
                :default_producer_indices;
                default=0,
            ),
        ),
        ParameterDefinition(
            ParameterSpec(
                :nutrient_half_saturation,
                :vector;
                doc="Nutrient half-saturation constant per plankton class.",
            ),
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27),
                :default_producer_indices;
                default=0,
            ),
        ),
        ParameterDefinition(
            ParameterSpec(
                :alpha,
                :vector;
                doc="Initial slope of the P-I curve per plankton class.",
            ),
            DiameterIndexedVectorDefault(0.1953 / 86400, :default_producer_indices; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :maximum_predation_rate,
                :vector;
                doc="Maximum zooplankton grazing rate per plankton class.",
            ),
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16),
                :default_consumer_indices;
                default=0,
            ),
        ),
        ParameterDefinition(
            ParameterSpec(
                :holling_half_saturation,
                :vector;
                doc="Holling type II half-saturation constant per plankton class.",
            ),
            DiameterIndexedVectorDefault(5.0, :default_consumer_indices; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :palatability_matrix,
                :matrix;
                axes=(:consumer, :prey),
                doc="Preference of each consumer for each prey class.",
            ),
            NoDefault(),
        ),
        ParameterDefinition(
            ParameterSpec(
                :assimilation_matrix,
                :matrix;
                axes=(:consumer, :prey),
                doc="Assimilation efficiency of each consumer on each prey class.",
            ),
            NoDefault(),
        ),
        ParameterDefinition(
            ParameterSpec(
                :optimum_predator_prey_ratio,
                :vector;
                doc="Preferred predator:prey diameter ratio per consumer (used to derive palatability_matrix).",
            ),
            DiameterIndexedVectorDefault(10.0, :default_consumer_indices; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :specificity,
                :vector;
                doc="Unimodal palatability specificity per consumer (used to derive palatability_matrix).",
            ),
            DiameterIndexedVectorDefault(0.3, :default_consumer_indices; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :protection,
                :vector;
                doc="Prey protection factor (used to derive palatability_matrix).",
            ),
            DiameterIndexedVectorDefault(1.0, :default_consumer_indices; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :assimilation_efficiency,
                :vector;
                doc="Assimilation efficiency per consumer (used to derive assimilation_matrix).",
            ),
            DiameterIndexedVectorDefault(0.32, :default_consumer_indices; default=0),
        ),
    )
end

function matrix_definitions(::NiPiZDFactory)
    return (;
        palatability_matrix=MatrixDefinition(PalatabilityAllometric()),
        assimilation_matrix=MatrixDefinition(AssimilationBinary()),
    )
end
