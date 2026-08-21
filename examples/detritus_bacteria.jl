# # [Model with bacterioplankton] (@id detritus_bacteria)
#
# Material pools and living prey are separate concepts. Here a structured POM
# pool is consumed by heterotrophic bacterioplankton `B`, while zooplankton `Z`
# graze the living bacterial population. A Q10 factor modifies POM consumption.

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Parameters:
    ConstDefault, DerivedDefault, FillDefault, ParameterDefinition, ParameterProvision
using Agate.Introspection: auxiliary_field_names, tracer_names
using Agate.Processes:
    Consumption, Grazing, ModelDefinition, Temperature, normalize_model, participants
using Oceananigans.Units: day

nothing #hide

# ## Components
#
# `POM` expands into two material tracers. `B` and `Z` are ordinary living
# populations; their ecological roles are supplied by the processes below.

components = (
    N=Pool(:nitrogen),
    POM=Pool(:nitrogen; size_structure=[0.5, 5.0]),
    B=Population(; currency=:nitrogen, size_structure=[0.8]),
    Z=Population(; currency=:nitrogen, size_structure=[10.0]),
)

processes = (
    consume_POM=Consumption(
        :heterotrophic;
        consumer=:B,
        resource=:POM,
        unassimilated_destination=:N,
        factors=(temperature=Temperature(:q10),),
    ),
    graze_bacteria=Grazing(
        :preferential; consumer=:Z, resource=:B, unassimilated_destination=:N
    ),
)

# ## Parameters
#
# Living-class parameters use the explicit global `:plankton` storage axis.
# POM half-saturation and bacterial assimilation omit `axes`, so their storage
# follows the process-local POM and B-by-POM applicability directly.

parameters = (
    ParameterDefinition(
        :maximum_consumption_rate,
        FillDefault(0.8 / day);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:consume_POM, :maximum_rate),
    ),
    ParameterDefinition(
        :pom_half_saturation,
        FillDefault(0.2);
        shape=:vector,
        provides=ParameterProvision(:consume_POM, :half_saturation),
    ),
    ParameterDefinition(
        :bacterial_assimilation,
        FillDefault(0.65);
        shape=:matrix,
        provides=ParameterProvision(:consume_POM, :assimilation),
    ),
    ParameterDefinition(
        :temperature_q10,
        ConstDefault(2.0);
        provides=ParameterProvision(:consume_POM, :q10),
    ),
    ParameterDefinition(
        :reference_temperature,
        ConstDefault(20.0);
        provides=ParameterProvision(:consume_POM, :reference_temperature),
    ),
    ParameterDefinition(
        :maximum_predation_rate,
        FillDefault(0.6 / day);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:graze_bacteria, :maximum_rate),
    ),
    ParameterDefinition(
        :holling_half_saturation,
        FillDefault(0.1);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:graze_bacteria, :half_saturation),
    ),
    ParameterDefinition(
        :optimum_predator_prey_ratio,
        FillDefault(12.5);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:graze_bacteria, :optimum_predator_prey_ratio),
    ),
    ParameterDefinition(
        :specificity,
        FillDefault(0.4);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:graze_bacteria, :specificity),
    ),
    ParameterDefinition(
        :protection,
        FillDefault(0.0);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:graze_bacteria, :protection),
    ),
    ParameterDefinition(
        :assimilation_efficiency,
        FillDefault(0.7);
        shape=:vector,
        axes=:plankton,
        provides=ParameterProvision(:graze_bacteria, :assimilation_efficiency),
    ),
    ParameterDefinition(
        :living_palatability,
        DerivedDefault(
            PalatabilityAllometric();
            deps=(:optimum_predator_prey_ratio, :specificity, :protection),
        );
        shape=:matrix,
        axes=(:consumer, :prey),
        runtime_path=(:interactions, :living_palatability),
        provides=ParameterProvision(:graze_bacteria, :palatability),
    ),
    ParameterDefinition(
        :living_assimilation,
        DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,));
        shape=:matrix,
        axes=(:consumer, :prey),
        runtime_path=(:interactions, :living_assimilation),
        provides=ParameterProvision(:graze_bacteria, :assimilation),
    ),
)

# ## Construct and inspect

definition = ModelDefinition(; components, processes, parameters)
normalized = normalize_model(definition)
bgc = construct(definition)

println("tracers: ", tracer_names(bgc))
println("drivers: ", auxiliary_field_names(bgc))
println("POM consumer: ", participants(normalized.processes.consume_POM).consumer)
println("living grazing prey: ", participants(normalized.processes.graze_bacteria).resource)

nothing #hide
