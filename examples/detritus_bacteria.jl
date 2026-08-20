# # [Detritus and bacterioplankton] (@id detritus_bacteria)
#
# Material pools and living prey are separate concepts. Here a structured POM
# pool is consumed by heterotrophic bacterioplankton `B`, while zooplankton `Z`
# graze the living bacterial population. A Q10 factor modifies POM consumption.

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Factories:
    ConstDefault, DerivedDefault, FillDefault,
    ParameterDefinition, ParameterProvision, ParameterSpec
using Agate.Introspection: auxiliary_field_names, tracer_names
using Agate.Processes:
    Consumption, Grazing, ModelDefinition, Temperature, normalize_model, participants
using Oceananigans.Units: day

nothing #hide

provision(process, path, formulation, slot; qualifier=NamedTuple()) =
    ParameterProvision(process, path, formulation, slot; qualifier)

function parameter(name, shape, provider, provides; axes=nothing)
    return ParameterDefinition(ParameterSpec(name, shape; axes, provides), provider)
end

scalar_parameter(name, value, provides) =
    parameter(name, :scalar, ConstDefault(value), provides)

vector_parameter(name, value, provides) =
    parameter(name, :vector, FillDefault(value), provides; axes=:plankton)

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
        :preferential;
        consumer=:Z,
        resource=:B,
        unassimilated_destination=:N,
    ),
)

# ## Parameters

parameters = (
    vector_parameter(
        :maximum_consumption_rate,
        0.8 / day,
        provision(:consume_POM, (), :heterotrophic, :maximum_rate),
    ),
    parameter(
        :pom_half_saturation,
        :vector,
        FillDefault(0.2),
        provision(:consume_POM, (), :heterotrophic, :half_saturation),
    ),
    parameter(
        :bacterial_assimilation,
        :matrix,
        FillDefault(0.65),
        provision(:consume_POM, (), :heterotrophic, :assimilation),
    ),
    scalar_parameter(
        :temperature_q10,
        2.0,
        provision(:consume_POM, (:factors, :temperature), :q10, :q10),
    ),
    scalar_parameter(
        :reference_temperature,
        20.0,
        provision(
            :consume_POM,
            (:factors, :temperature),
            :q10,
            :reference_temperature,
        ),
    ),
    vector_parameter(
        :maximum_predation_rate,
        0.6 / day,
        provision(:graze_bacteria, (), :preferential, :maximum_rate),
    ),
    vector_parameter(
        :holling_half_saturation,
        0.1,
        provision(:graze_bacteria, (), :preferential, :half_saturation),
    ),
    vector_parameter(
        :optimum_predator_prey_ratio,
        12.5,
        provision(
            :graze_bacteria,
            (:palatability, :default),
            :allometric,
            :optimum_predator_prey_ratio,
        ),
    ),
    vector_parameter(
        :specificity,
        0.4,
        provision(
            :graze_bacteria,
            (:palatability, :default),
            :allometric,
            :specificity,
        ),
    ),
    vector_parameter(
        :protection,
        0.0,
        provision(
            :graze_bacteria,
            (:palatability, :default),
            :allometric,
            :protection,
        ),
    ),
    vector_parameter(
        :assimilation_efficiency,
        0.7,
        provision(
            :graze_bacteria,
            (:assimilation, :default),
            :binary,
            :assimilation_efficiency,
        ),
    ),
    parameter(
        :living_palatability,
        :matrix,
        DerivedDefault(
            PalatabilityAllometric();
            deps=(:optimum_predator_prey_ratio, :specificity, :protection),
        ),
        provision(:graze_bacteria, (), :preferential, :palatability);
        axes=(:consumer, :prey),
    ),
    parameter(
        :living_assimilation,
        :matrix,
        DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,)),
        provision(:graze_bacteria, (), :preferential, :assimilation);
        axes=(:consumer, :prey),
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
