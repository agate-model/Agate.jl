# # [Model with bacterioplankton] (@id detritus_bacteria)
#
# Material pools and living prey are separate concepts. Here a structured POM
# pool is consumed by heterotrophic bacterioplankton `B`, while zooplankton `Z`
# graze the living bacterial population. A Q10 factor modifies POM consumption.

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Parameters:
    ConstantDefault, DerivedDefault, Parameter, ParameterProvision
using Agate.Introspection: auxiliary_field_names, tracer_names
using Agate.Processes:
    Consumption, ModelDefinition, ProductRouting, Temperature, Q10,
    HeterotrophicConsumption, PreferentialGrazing, DirectRouting, normalize_model, participants
using Oceananigans.Units: day

nothing #hide

# ## Components
#
# `POM` expands into two material tracers. `B` and `Z` are ordinary living
# populations; their ecological roles are supplied by the processes below.

components = (
    N=Pool(:nitrogen),
    POM=Pool(:nitrogen; size_structure=[0.5, 5.0]),
    B=Population(:nitrogen; size_structure=[0.8]),
    Z=Population(:nitrogen; size_structure=[10.0]),
)

processes = (
    consume_POM=Consumption(
        HeterotrophicConsumption();
        consumers=:B,
        resources=:POM,
        factors=(temperature=Temperature(Q10()),),
        routing=ProductRouting(DirectRouting(); destination=:N),
    ),
    graze_bacteria=Consumption(
        PreferentialGrazing();
        consumers=:Z,
        resources=:B,
        routing=ProductRouting(DirectRouting(); destination=:N),
    ),
)

# ## Parameters
#
# Each NamedTuple key is the stable model parameter name. Living-class parameters
# use the explicit global `:plankton` storage axis.
# POM half-saturation and bacterial assimilation omit `axes`, so their storage
# follows the process-local POM and B-by-POM applicability directly.

parameters = (
    maximum_consumption_rate=Parameter(
        ConstantDefault(0.8 / day);
        axes=:plankton,
        provides=ParameterProvision(:consume_POM, :maximum_rate),
    ),
    pom_half_saturation=Parameter(
        ConstantDefault(0.2);
        provides=ParameterProvision(:consume_POM, :half_saturation),
    ),
    bacterial_assimilation=Parameter(
        ConstantDefault(0.65);
        provides=ParameterProvision(:consume_POM, :assimilation),
    ),
    temperature_q10=Parameter(
        ConstantDefault(2.0);
        provides=ParameterProvision(:consume_POM, :q10),
    ),
    reference_temperature=Parameter(
        ConstantDefault(20.0);
        provides=ParameterProvision(:consume_POM, :reference_temperature),
    ),
    maximum_predation_rate=Parameter(
        ConstantDefault(0.6 / day);
        axes=:plankton,
        provides=ParameterProvision(:graze_bacteria, :maximum_rate),
    ),
    holling_half_saturation=Parameter(
        ConstantDefault(0.1);
        axes=:plankton,
        provides=ParameterProvision(:graze_bacteria, :half_saturation),
    ),
    optimum_predator_prey_ratio=Parameter(
        ConstantDefault(12.5);
        axes=:plankton,
    ),
    specificity=Parameter(
        ConstantDefault(0.4);
        axes=:plankton,
    ),
    protection=Parameter(
        ConstantDefault(0.0);
        axes=:plankton,
    ),
    assimilation_efficiency=Parameter(
        ConstantDefault(0.7);
        axes=:plankton,
    ),
    living_palatability=Parameter(
        DerivedDefault(
            PalatabilityAllometric();
            deps=(:optimum_predator_prey_ratio, :specificity, :protection),
        );
        axes=(:consumer, :prey),
        provides=ParameterProvision(:graze_bacteria, :palatability),
    ),
    living_assimilation=Parameter(
        DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,));
        axes=(:consumer, :prey),
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
