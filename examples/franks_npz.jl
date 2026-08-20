# # [Franks et al. (1986) NPZ] (@id franks_npz)
#
# Agate models are defined from structured components, named biological processes,
# and parameter definitions. These pieces form a `ModelDefinition`; `construct`
# realizes concrete tracers, resolves parameters and drivers, and compiles the
# runtime equations. This NPZ example is the smallest complete introduction to
# that pattern.
#
# A compact `ModelDefinition` is enough to describe a nutrient--phytoplankton--
# zooplankton ecosystem. The classic model of
# [Franks, Wroblewski & Flierl (1986)](https://doi.org/10.1007/BF00397577)
# has three nitrogen state variables, Michaelis--Menten nutrient uptake,
# herbivore grazing, mortality, and immediate recycling to dissolved nutrient.
#
# The source paper compares Ivlev and food-acclimating herbivore responses.
# This tutorial keeps the three-state NPZ bookkeeping and translates grazing to
# Agate's preferential (Holling-II) formulation.

using Agate.Configuration: AssimilationBinary, PalatabilityAllometric, Population, Pool
using Agate.Construction: construct
using Agate.Factories:
    ConstDefault, DerivedDefault, FillDefault,
    ParameterDefinition, ParameterProvision, ParameterSpec
using Agate.Introspection: auxiliary_field_names, tracer_names
using Agate.Processes:
    Grazing, Growth, Light, ModelDefinition, Mortality,
    NutrientResponse, ProductRouting
using Oceananigans.Units: day

nothing #hide

# ## Parameter declarations
#
# Process formulations declare scientific requirements such as a maximum growth
# rate or grazing assimilation. `ParameterProvision` binds each requirement to a
# named model parameter.

provision(process, path, formulation, slot; qualifier=NamedTuple()) =
    ParameterProvision(process, path, formulation, slot; qualifier)

function parameter(name, shape, provider, provides; axes=nothing)
    spec = ParameterSpec(name, shape; axes, provides)
    return ParameterDefinition(spec, provider)
end

scalar_parameter(name, value, provides) =
    parameter(name, :scalar, ConstDefault(value), provides)

vector_parameter(name, value, provides) =
    parameter(name, :vector, FillDefault(value), provides; axes=:plankton)

nothing #hide

# ## Components and processes
#
# Components carry structural identity. Ecological roles come from participation
# in named processes.

components = (
    N=Pool(:nitrogen),
    P=Population(; currency=:nitrogen, size_structure=[1.0]),
    Z=Population(; currency=:nitrogen, size_structure=[10.0]),
)

recycle_to_N = ProductRouting(:partition; retained=:N, exported=:N)

processes = (
    growth_P=Growth(;
        population=:P,
        factors=(
            light=Light(:smith; driver=:PAR),
            nutrients=NutrientResponse(:monod; resource=:N),
        ),
    ),
    grazing_Z_on_P=Grazing(
        :preferential;
        consumer=:Z,
        resource=:P,
        unassimilated_destination=:N,
    ),
    mortality_P=Mortality(:linear; population=:P, routing=recycle_to_N),
    mortality_Z=Mortality(:linear; population=:Z, routing=recycle_to_N),
)

# The numerical defaults below provide a compact teaching parameterization.
# A strong Smith slope makes `PAR` effectively non-limiting when a constant
# order-one light driver is supplied. The Holling half-saturation belongs to
# Agate's preferential grazing formulation.

parameters = (
    vector_parameter(
        :maximum_growth_rate,
        0.69 / day,
        provision(:growth_P, (:factors, :light), :smith, :maximum_rate),
    ),
    vector_parameter(
        :alpha,
        100 / day,
        provision(:growth_P, (:factors, :light), :smith, :alpha),
    ),
    vector_parameter(
        :nutrient_half_saturation,
        1.0,
        provision(
            :growth_P,
            (:factors, :nutrients),
            :monod,
            :K;
            qualifier=(resource=:N,),
        ),
    ),
    vector_parameter(
        :maximum_predation_rate,
        1.5 / day,
        provision(:grazing_Z_on_P, (), :preferential, :maximum_rate),
    ),
    vector_parameter(
        :holling_half_saturation,
        1.0,
        provision(:grazing_Z_on_P, (), :preferential, :half_saturation),
    ),
    vector_parameter(
        :optimum_predator_prey_ratio,
        10.0,
        provision(
            :grazing_Z_on_P,
            (:palatability, :default),
            :allometric,
            :optimum_predator_prey_ratio,
        ),
    ),
    vector_parameter(
        :specificity,
        0.5,
        provision(
            :grazing_Z_on_P,
            (:palatability, :default),
            :allometric,
            :specificity,
        ),
    ),
    vector_parameter(
        :protection,
        0.0,
        provision(
            :grazing_Z_on_P,
            (:palatability, :default),
            :allometric,
            :protection,
        ),
    ),
    vector_parameter(
        :assimilation_efficiency,
        0.7,
        provision(
            :grazing_Z_on_P,
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
        provision(:grazing_Z_on_P, (), :preferential, :palatability);
        axes=(:consumer, :prey),
    ),
    parameter(
        :assimilation_matrix,
        :matrix,
        DerivedDefault(AssimilationBinary(); deps=(:assimilation_efficiency,)),
        provision(:grazing_Z_on_P, (), :preferential, :assimilation);
        axes=(:consumer, :prey),
    ),
    vector_parameter(
        :phytoplankton_mortality,
        0.1 / day,
        provision(
            :mortality_P,
            (),
            :linear,
            :rate;
            qualifier=(population=:P,),
        ),
    ),
    vector_parameter(
        :zooplankton_mortality,
        0.2 / day,
        provision(
            :mortality_Z,
            (),
            :linear,
            :rate;
            qualifier=(population=:Z,),
        ),
    ),
    scalar_parameter(
        :phytoplankton_mortality_export_fraction,
        0.0,
        provision(:mortality_P, (:routing,), :partition, :export_fraction),
    ),
    scalar_parameter(
        :zooplankton_mortality_export_fraction,
        0.0,
        provision(:mortality_Z, (:routing,), :partition, :export_fraction),
    ),
)

# ## Construct directly

definition = ModelDefinition(; components, processes, parameters)
bgc = construct(definition)

println("tracers: ", tracer_names(bgc))
println("drivers: ", auxiliary_field_names(bgc))
println("maximum growth: ", bgc.parameters.maximum_growth_rate[1] * day, " day^-1")

nothing #hide
