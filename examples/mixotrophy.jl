# # [Model with mixotrophy] (@id mixotrophy)
#
# A mixotroph does not need a special component type. `M` is an ordinary
# `Plankton` that participates in both growth and grazing.

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Plankton, Pool
using Agate.Construction: construct
using Agate.Parameters: DerivedDefault, MetaParameter, Parameter
using Agate.Introspection: auxiliary_field_names, tracer_names
using Agate.Processes:
    Consumption, Growth, Light, ModelDefinition, NutrientResponse, Smith, Monod, PreferentialGrazing, participants
using Oceananigans.Units: day

nothing #hide

components = (
    N=Pool(:nitrogen),
    P=Plankton(; states=:nitrogen, reference_state=:nitrogen, size_structure=[1.0]),
    M=Plankton(; states=:nitrogen, reference_state=:nitrogen, size_structure=[4.0]),
)

processes = (
    growth_M=Growth(;
        plankton=:M,
        bindings=(maximum_rate=:maximum_growth_rate,),
        factors=(
            light=Light(Smith(); driver=:PAR),
            nutrients=NutrientResponse(
                Monod(); resource=:N, bindings=(half_saturation=:nutrient_half_saturation,)
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
# their local scientific slots directly to these keys, and those slots determine runtime
# storage automatically. The interaction traits retain an explicit `:plankton` axis because
# they are dependency-only inputs to the derived interaction matrices.

parameters = (
    maximum_growth_rate=Parameter(0.8 / day),
    alpha=Parameter(0.08 / day),
    nutrient_half_saturation=Parameter(0.2),
    maximum_predation_rate=Parameter(0.4 / day),
    holling_half_saturation=Parameter(0.15),
    optimum_predator_prey_ratio=MetaParameter(4.0; axes=:plankton),
    specificity=MetaParameter(0.5; axes=:plankton),
    protection=MetaParameter(0.0; axes=:plankton),
    assimilation_efficiency=MetaParameter(0.65; axes=:plankton),
    palatability_matrix=Parameter(
        DerivedDefault(
            PalatabilityAllometric();
            deps=(:optimum_predator_prey_ratio, :specificity, :protection),
        )
    ),
    assimilation_matrix=Parameter(
        DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,))
    ),
)

definition = ModelDefinition(; components, processes, parameters)
bgc = construct(definition)

println("tracers: ", tracer_names(bgc))
println("drivers: ", auxiliary_field_names(bgc))
println("M grows as: ", participants(processes.growth_M).plankton)
println("M grazes as: ", participants(processes.grazing_M).consumer)

nothing #hide
