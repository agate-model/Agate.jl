# # [Model with mixotrophy] (@id mixotrophy)
#
# A mixotroph does not need a special component type. `M` is an ordinary
# `Population` that participates in both growth and grazing.

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Parameters:
    DerivedDefault, FillDefault, ParameterDefinition, ParameterProvision
using Agate.Introspection: auxiliary_field_names, tracer_names
using Agate.Processes:
    Grazing, Growth, Light, ModelDefinition, NutrientResponse, normalize_model, participants
using Oceananigans.Units: day

nothing #hide

components = (
    N=Pool(:nitrogen),
    P=Population(; currency=:nitrogen, size_structure=[1.0]),
    M=Population(; currency=:nitrogen, size_structure=[4.0]),
)

processes = (
    growth_M=Growth(;
        population=:M,
        factors=(
            light=Light(:smith; driver=:PAR),
            nutrients=NutrientResponse(:monod; resource=:N),
        ),
    ),
    grazing_M=Grazing(
        :preferential; consumer=:M, resource=:P, unassimilated_destination=:N
    ),
)

# ## Parameters
#
# Parameters name their scientific process slot through `ParameterProvision`;
# formulation and nested path are inferred unless disambiguation is needed.
# Vector parameters use the realized plankton axis, while derived interaction
# matrices are stored under `bgc.parameters.interactions`.

parameters = (
    ParameterDefinition(
        :maximum_growth_rate,
        FillDefault(0.8 / day);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:growth_M, :maximum_rate),
    ),
    ParameterDefinition(
        :alpha,
        FillDefault(0.08 / day);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:growth_M, :alpha),
    ),
    ParameterDefinition(
        :nutrient_half_saturation,
        FillDefault(0.2);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:growth_M, :K),
    ),
    ParameterDefinition(
        :maximum_predation_rate,
        FillDefault(0.4 / day);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:grazing_M, :maximum_rate),
    ),
    ParameterDefinition(
        :holling_half_saturation,
        FillDefault(0.15);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:grazing_M, :half_saturation),
    ),
    ParameterDefinition(
        :optimum_predator_prey_ratio,
        FillDefault(4.0);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:grazing_M, :optimum_predator_prey_ratio),
    ),
    ParameterDefinition(
        :specificity,
        FillDefault(0.5);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:grazing_M, :specificity),
    ),
    ParameterDefinition(
        :protection,
        FillDefault(0.0);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:grazing_M, :protection),
    ),
    ParameterDefinition(
        :assimilation_efficiency,
        FillDefault(0.65);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:grazing_M, :assimilation_efficiency),
    ),
    ParameterDefinition(
        :palatability_matrix,
        DerivedDefault(
            PalatabilityAllometric();
            deps=(:optimum_predator_prey_ratio, :specificity, :protection),
        );
        shape=:matrix,
        axes=(:consumer, :prey),
        runtime_path=(:interactions, :palatability),
        provides=ParameterProvision(:grazing_M, :palatability),
    ),
    ParameterDefinition(
        :assimilation_matrix,
        DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,));
        shape=:matrix,
        axes=(:consumer, :prey),
        runtime_path=(:interactions, :assimilation),
        provides=ParameterProvision(:grazing_M, :assimilation),
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
