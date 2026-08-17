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
    DiameterIndexedVectorDefault,
    DiameterIndexedMaterialization

import ...Configuration: matrix_definitions

using ...Library.Allometry: AllometricParam, PowerLaw

using ...Configuration: MatrixDefinition, PalatabilityAllometric, AssimilationBinary

function parameter_definitions(::NiPiZDFactory)
    detritus_remin = 0.1213 / 86400
    all_plankton_materialization = DiameterIndexedMaterialization(; fill_value=0)
    producer_materialization =
        DiameterIndexedMaterialization(:producers; fill_value=0)
    consumer_materialization =
        DiameterIndexedMaterialization(:consumers; fill_value=0)

    return (
        ParameterDefinition(
            ParameterSpec(
                :detritus_remineralization, :scalar; doc="Detritus remineralization rate."
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
                axes=:plankton,
                materialization=all_plankton_materialization,
                doc="Linear mortality coefficient per plankton class.",
            ),
            FillDefault(8e-7),
        ),
        ParameterDefinition(
            ParameterSpec(
                :quadratic_mortality,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                doc="Quadratic mortality coefficient per plankton class.",
            ),
            DiameterIndexedVectorDefault(1e-6, :consumers; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :maximum_growth_rate,
                :vector;
                axes=:plankton,
                materialization=producer_materialization,
                doc="Maximum phytoplankton growth rate per plankton class.",
            ),
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15),
                :producers;
                default=0,
            ),
        ),
        ParameterDefinition(
            ParameterSpec(
                :nutrient_half_saturation,
                :vector;
                axes=:plankton,
                materialization=producer_materialization,
                doc="Nutrient half-saturation constant per plankton class.",
            ),
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27),
                :producers;
                default=0,
            ),
        ),
        ParameterDefinition(
            ParameterSpec(
                :alpha,
                :vector;
                axes=:plankton,
                materialization=producer_materialization,
                doc="Initial slope of the P-I curve per plankton class.",
            ),
            DiameterIndexedVectorDefault(
                0.1953 / 86400, :producers; default=0
            ),
        ),
        ParameterDefinition(
            ParameterSpec(
                :maximum_predation_rate,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                doc="Maximum zooplankton grazing rate per plankton class.",
            ),
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16),
                :consumers;
                default=0,
            ),
        ),
        ParameterDefinition(
            ParameterSpec(
                :holling_half_saturation,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                doc="Holling type II half-saturation constant per plankton class.",
            ),
            DiameterIndexedVectorDefault(5.0, :consumers; default=0),
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
                axes=:plankton,
                materialization=consumer_materialization,
                doc="Preferred predator:prey diameter ratio per consumer (used to derive palatability_matrix).",
            ),
            DiameterIndexedVectorDefault(10.0, :consumers; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :specificity,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                doc="Unimodal palatability specificity per consumer (used to derive palatability_matrix).",
            ),
            DiameterIndexedVectorDefault(0.3, :consumers; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :protection,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                doc="Prey protection factor (used to derive palatability_matrix).",
            ),
            DiameterIndexedVectorDefault(1.0, :consumers; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :assimilation_efficiency,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                doc="Assimilation efficiency per consumer (used to derive assimilation_matrix).",
            ),
            DiameterIndexedVectorDefault(0.32, :consumers; default=0),
        ),
    )
end

function matrix_definitions(::NiPiZDFactory)
    return (;
        palatability_matrix=MatrixDefinition(PalatabilityAllometric()),
        assimilation_matrix=MatrixDefinition(AssimilationBinary()),
    )
end
