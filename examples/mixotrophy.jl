# # [Model with mixotrophy] (@id mixotrophy)
#
# A mixotroph does not need a special component type. `M` is an ordinary
# `Population` that participates in both growth and grazing.

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Factories:
    DerivedDefault, FillDefault,
    ParameterDefinition, ParameterProvision, ParameterSpec
using Agate.Introspection: auxiliary_field_names, tracer_names
using Agate.Processes:
    Grazing, Growth, Light, ModelDefinition, NutrientResponse,
    normalize_model, participants
using Oceananigans.Units: day

nothing #hide

provision(process, path, formulation, slot; qualifier=NamedTuple()) =
    ParameterProvision(process, path, formulation, slot; qualifier)

function parameter(name, shape, provider, provides; axes=nothing)
    return ParameterDefinition(ParameterSpec(name, shape; axes, provides), provider)
end

vector_parameter(name, value, provides) =
    parameter(name, :vector, FillDefault(value), provides; axes=:plankton)

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
        :preferential;
        consumer=:M,
        resource=:P,
        unassimilated_destination=:N,
    ),
)

parameters = (
    vector_parameter(
        :maximum_growth_rate,
        0.8 / day,
        provision(:growth_M, (:factors, :light), :smith, :maximum_rate),
    ),
    vector_parameter(
        :alpha,
        0.08 / day,
        provision(:growth_M, (:factors, :light), :smith, :alpha),
    ),
    vector_parameter(
        :nutrient_half_saturation,
        0.2,
        provision(
            :growth_M,
            (:factors, :nutrients),
            :monod,
            :K;
            qualifier=(resource=:N,),
        ),
    ),
    vector_parameter(
        :maximum_predation_rate,
        0.4 / day,
        provision(:grazing_M, (), :preferential, :maximum_rate),
    ),
    vector_parameter(
        :holling_half_saturation,
        0.15,
        provision(:grazing_M, (), :preferential, :half_saturation),
    ),
    vector_parameter(
        :optimum_predator_prey_ratio,
        4.0,
        provision(
            :grazing_M,
            (:palatability, :default),
            :allometric,
            :optimum_predator_prey_ratio,
        ),
    ),
    vector_parameter(
        :specificity,
        0.5,
        provision(
            :grazing_M,
            (:palatability, :default),
            :allometric,
            :specificity,
        ),
    ),
    vector_parameter(
        :protection,
        0.0,
        provision(
            :grazing_M,
            (:palatability, :default),
            :allometric,
            :protection,
        ),
    ),
    vector_parameter(
        :assimilation_efficiency,
        0.65,
        provision(
            :grazing_M,
            (:assimilation, :default),
            :binary,
            :assimilation_efficiency,
        ),
    ),
    parameter(
        :palatability_matrix,
        :matrix,
        DerivedDefault(
            PalatabilityAllometric();
            deps=(:optimum_predator_prey_ratio, :specificity, :protection),
        ),
        provision(:grazing_M, (), :preferential, :palatability);
        axes=(:consumer, :prey),
    ),
    parameter(
        :assimilation_matrix,
        :matrix,
        DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,)),
        provision(:grazing_M, (), :preferential, :assimilation);
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
