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
    ParameterProvision,
    ConstDefault,
    NoDefault,
    FillDefault,
    DiameterIndexedVectorDefault,
    DiameterIndexedMaterialization

import ...Configuration: matrix_definitions

using ...Library.Allometry: AllometricParam, PowerLaw

using ...Configuration: MatrixDefinition, PalatabilityAllometric, AssimilationBinary

provision(process, path, formulation, slot; qualifier=NamedTuple()) =
    ParameterProvision(process, path, formulation, slot; qualifier)

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
                :detritus_remineralization,
                :scalar;
                provides=provision(
                    :remineralization_D,
                    (),
                    :linear,
                    :rate;
                    qualifier=(source=:D,),
                ),
            ),
            ConstDefault(detritus_remin),
        ),
        ParameterDefinition(
            ParameterSpec(
                :mortality_export_fraction,
                :scalar;
                provides=(
                    provision(
                        :linear_mortality_P, (:routing,), :partition, :export_fraction
                    ),
                    provision(
                        :linear_mortality_Z, (:routing,), :partition, :export_fraction
                    ),
                    provision(
                        :quadratic_mortality_Z, (:routing,), :partition, :export_fraction
                    ),
                ),
            ),
            ConstDefault(0.2),
        ),
        ParameterDefinition(
            ParameterSpec(
                :linear_mortality,
                :vector;
                axes=:plankton,
                materialization=all_plankton_materialization,
                provides=(
                    provision(
                        :linear_mortality_P,
                        (),
                        :linear,
                        :rate;
                        qualifier=(population=:P,),
                    ),
                    provision(
                        :linear_mortality_Z,
                        (),
                        :linear,
                        :rate;
                        qualifier=(population=:Z,),
                    ),
                ),
            ),
            FillDefault(8e-7),
        ),
        ParameterDefinition(
            ParameterSpec(
                :quadratic_mortality,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                provides=provision(
                    :quadratic_mortality_Z,
                    (),
                    :quadratic,
                    :rate;
                    qualifier=(population=:Z,),
                ),
            ),
            DiameterIndexedVectorDefault(1e-6, :consumers; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :maximum_growth_rate,
                :vector;
                axes=:plankton,
                materialization=producer_materialization,
                provides=provision(:growth_P, (), :smith, :maximum_rate),
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
                provides=provision(
                    :growth_P,
                    (:limitation,),
                    :monod,
                    :K;
                    qualifier=(resource=:N,),
                ),
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
                provides=provision(:growth_P, (), :smith, :alpha),
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
                provides=provision(
                    :grazing_Z_on_P, (), :preferential, :maximum_rate
                ),
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
                provides=provision(
                    :grazing_Z_on_P, (), :preferential, :half_saturation
                ),
            ),
            DiameterIndexedVectorDefault(5.0, :consumers; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :palatability_matrix,
                :matrix;
                axes=(:consumer, :prey),
                provides=provision(
                    :grazing_Z_on_P, (), :preferential, :palatability
                ),
            ),
            NoDefault(),
        ),
        ParameterDefinition(
            ParameterSpec(
                :assimilation_matrix,
                :matrix;
                axes=(:consumer, :prey),
                provides=provision(
                    :grazing_Z_on_P, (), :preferential, :assimilation
                ),
            ),
            NoDefault(),
        ),
        ParameterDefinition(
            ParameterSpec(
                :optimum_predator_prey_ratio,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                provides=provision(
                    :grazing_Z_on_P,
                    (:palatability, :default),
                    :allometric,
                    :optimum_predator_prey_ratio,
                ),
            ),
            DiameterIndexedVectorDefault(10.0, :consumers; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :specificity,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                provides=provision(
                    :grazing_Z_on_P, (:palatability, :default), :allometric, :specificity
                ),
            ),
            DiameterIndexedVectorDefault(0.3, :consumers; default=0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :protection,
                :vector;
                axes=:plankton,
                materialization=producer_materialization,
                provides=provision(
                    :grazing_Z_on_P, (:palatability, :default), :allometric, :protection
                ),
            ),
            DiameterIndexedVectorDefault(0.0, :producers; default=1.0),
        ),
        ParameterDefinition(
            ParameterSpec(
                :assimilation_efficiency,
                :vector;
                axes=:plankton,
                materialization=consumer_materialization,
                provides=provision(
                    :grazing_Z_on_P,
                    (:assimilation, :default),
                    :binary,
                    :assimilation_efficiency,
                ),
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
