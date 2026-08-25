# # [Model with mixotrophy] (@id mixotrophy)
#
# A mixotroph does not need a special component type. `M` is an ordinary
# `Population` that participates in both growth and grazing.

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Parameters: DerivedDefault, ConstantDefault, Parameter
using Agate.Introspection: auxiliary_field_names, tracer_names
using Agate.Processes:
    Consumption, Growth, Light, ModelDefinition, NutrientResponse, Smith, Monod, PreferentialGrazing, normalize_model, participants
using Oceananigans.Units: day

nothing #hide

components = (
    N=Pool(:nitrogen),
    P=Population(:nitrogen; size_structure=[1.0]),
    M=Population(:nitrogen; size_structure=[4.0]),
)

processes = (
    growth_M=Growth(;
        populations=:M,
        factors=(
            light=Light(
                Smith(); driver=:PAR, bindings=(maximum_rate=:maximum_growth_rate,)
            ),
            nutrients=NutrientResponse(
                Monod(); resource=:N, bindings=(K=:nutrient_half_saturation,)
            ),
        ),
    ),
    grazing_M=Consumption(
        PreferentialGrazing();
        consumers=:M,
        resources=:P,
        bindings=(
            maximum_rate=:maximum_predation_rate,
            half_saturation=:holling_half_saturation,
            palatability=:palatability_matrix,
            assimilation=:assimilation_matrix,
        ),
        unassimilated_products=:N,
    ),
)

# ## Parameters
#
# Each NamedTuple key is the stable model parameter name. Processes and factors bind
# their local scientific slots directly to these keys. Vector parameters use the
# realized plankton axis, while derived interaction matrices are ordinary top-level
# runtime parameters.

parameters = (
    maximum_growth_rate=Parameter(
        ConstantDefault(0.8 / day);
        axes=:plankton,
    ),
    alpha=Parameter(
        ConstantDefault(0.08 / day);
        axes=:plankton,
    ),
    nutrient_half_saturation=Parameter(
        ConstantDefault(0.2);
        axes=:plankton,
    ),
    maximum_predation_rate=Parameter(
        ConstantDefault(0.4 / day);
        axes=:plankton,
    ),
    holling_half_saturation=Parameter(
        ConstantDefault(0.15);
        axes=:plankton,
    ),
    optimum_predator_prey_ratio=Parameter(
        ConstantDefault(4.0);
        axes=:plankton,
    ),
    specificity=Parameter(
        ConstantDefault(0.5);
        axes=:plankton,
    ),
    protection=Parameter(
        ConstantDefault(0.0);
        axes=:plankton,
    ),
    assimilation_efficiency=Parameter(
        ConstantDefault(0.65);
        axes=:plankton,
    ),
    palatability_matrix=Parameter(
        DerivedDefault(
            PalatabilityAllometric();
            deps=(:optimum_predator_prey_ratio, :specificity, :protection),
        );
        axes=(:consumer, :prey),
    ),
    assimilation_matrix=Parameter(
        DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,));
        axes=(:consumer, :prey),
    ),
)

definition = ModelDefinition(; components, processes, parameters)
normalized = normalize_model(definition)
bgc = construct(definition)

println("tracers: ", tracer_names(bgc))
println("drivers: ", auxiliary_field_names(bgc))
println("M grows as: ", participants(normalized.processes.growth_M).population)
println("M grazes as: ", participants(normalized.processes.grazing_M).consumer)

nothing #hide
