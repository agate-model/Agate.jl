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
    ConstDefault,
    DerivedDefault,
    FillDefault,
    DiameterIndexedVectorDefault,
    DiameterIndexedMaterialization

using ...Library.Allometry: AllometricParam, PowerLaw

using ...Configuration: PalatabilityAllometric, AssimilationBinary

function parameter_definitions(::NiPiZDFamily)
    detritus_remin = 0.1213 / 86400
    all_plankton_materialization = DiameterIndexedMaterialization(; fill_value=0)
    producer_materialization =
        DiameterIndexedMaterialization(; fill_value=0)
    consumer_materialization =
        DiameterIndexedMaterialization(; fill_value=0)

    return (
        ParameterDefinition(
            :detritus_remineralization,
            ConstDefault(detritus_remin);
            provides=ParameterProvision(:remineralization_D, :linear, :rate; qualifier=(source=:D,),),
        ),
        ParameterDefinition(
            :mortality_export_fraction,
            ConstDefault(0.2);
            provides=(
                ParameterProvision(:linear_mortality_P, :partition, :export_fraction; path=(:routing,)),
                ParameterProvision(:linear_mortality_Z, :partition, :export_fraction; path=(:routing,)),
                ParameterProvision(:quadratic_mortality_Z, :partition, :export_fraction; path=(:routing,)),
            ),
        ),
        ParameterDefinition(
            :linear_mortality,
            FillDefault(8e-7);
            shape=:vector,
            axes=:plankton,
            materialization=all_plankton_materialization,
            provides=(
                ParameterProvision(:linear_mortality_P, :linear, :rate; qualifier=(population=:P,),),
                ParameterProvision(:linear_mortality_Z, :linear, :rate; qualifier=(population=:Z,),),
            ),
        ),
        ParameterDefinition(
            :quadratic_mortality,
            DiameterIndexedVectorDefault(1e-6; default=0);
            shape=:vector,
            axes=:plankton,
            materialization=consumer_materialization,
            provides=ParameterProvision(:quadratic_mortality_Z, :quadratic, :rate; qualifier=(population=:Z,),),
        ),
        ParameterDefinition(
            :maximum_growth_rate,
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=2 / 86400, exponent=-0.15);
                default=0,
            );
            shape=:vector,
            axes=:plankton,
            materialization=producer_materialization,
            provides=ParameterProvision(:growth_P, :smith, :maximum_rate; path=(:factors, :light)),
        ),
        ParameterDefinition(
            :nutrient_half_saturation,
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=0.17, exponent=0.27);
                default=0,
            );
            shape=:vector,
            axes=:plankton,
            materialization=producer_materialization,
            provides=ParameterProvision(:growth_P, :monod, :K; path=(:factors, :nutrients), qualifier=(resource=:N,),),
        ),
        ParameterDefinition(
            :alpha,
            DiameterIndexedVectorDefault(
                0.1953 / 86400; default=0
            );
            shape=:vector,
            axes=:plankton,
            materialization=producer_materialization,
            provides=ParameterProvision(:growth_P, :smith, :alpha; path=(:factors, :light)),
        ),
        ParameterDefinition(
            :maximum_predation_rate,
            DiameterIndexedVectorDefault(
                AllometricParam(PowerLaw(); prefactor=30.84 / 86400, exponent=-0.16);
                default=0,
            );
            shape=:vector,
            axes=:plankton,
            materialization=consumer_materialization,
            provides=ParameterProvision(:grazing_Z_on_P, :preferential, :maximum_rate),
        ),
        ParameterDefinition(
            :holling_half_saturation,
            DiameterIndexedVectorDefault(5.0; default=0);
            shape=:vector,
            axes=:plankton,
            materialization=consumer_materialization,
            provides=ParameterProvision(:grazing_Z_on_P, :preferential, :half_saturation),
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
            shape=:matrix,
            axes=(:consumer, :prey),
            runtime_path=(:interactions, :palatability),
            provides=ParameterProvision(:grazing_Z_on_P, :preferential, :palatability),
        ),
        ParameterDefinition(
            :assimilation_matrix,
            DerivedDefault(
                AssimilationBinary(); deps=(:assimilation_efficiency,)
            );
            shape=:matrix,
            axes=(:consumer, :prey),
            runtime_path=(:interactions, :assimilation),
            provides=ParameterProvision(:grazing_Z_on_P, :preferential, :assimilation),
        ),
        ParameterDefinition(
            :optimum_predator_prey_ratio,
            DiameterIndexedVectorDefault(10.0; default=0);
            shape=:vector,
            axes=:plankton,
            materialization=consumer_materialization,
            provides=ParameterProvision(:grazing_Z_on_P, :allometric, :optimum_predator_prey_ratio; path=(:palatability, :default)),
        ),
        ParameterDefinition(
            :specificity,
            DiameterIndexedVectorDefault(0.3; default=0);
            shape=:vector,
            axes=:plankton,
            materialization=consumer_materialization,
            provides=ParameterProvision(:grazing_Z_on_P, :allometric, :specificity; path=(:palatability, :default)),
        ),
        ParameterDefinition(
            :protection,
            DiameterIndexedVectorDefault(0.0; default=1.0);
            shape=:vector,
            axes=:plankton,
            materialization=producer_materialization,
            provides=ParameterProvision(:grazing_Z_on_P, :allometric, :protection; path=(:palatability, :default)),
        ),
        ParameterDefinition(
            :assimilation_efficiency,
            DiameterIndexedVectorDefault(0.32; default=0);
            shape=:vector,
            axes=:plankton,
            materialization=consumer_materialization,
            provides=ParameterProvision(:grazing_Z_on_P, :binary, :assimilation_efficiency; path=(:assimilation, :default)),
        ),
    )
end
