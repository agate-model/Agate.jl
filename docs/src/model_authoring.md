# Model authoring

Agate models are authored from three pieces of scientific information:

1. **components** describe state identity and structure with [`Population`](@ref) and [`Pool`](@ref);
2. **processes** describe named biological transformations such as [`Growth`](@ref), [`Grazing`](@ref), [`Consumption`](@ref), and [`Mortality`](@ref);
3. **parameter definitions** bind named defaults to the semantic requirements declared by those process formulations.

These pieces form a [`ModelDefinition`](@ref), which can be constructed directly with [`Agate.Construction.construct`](@ref). Agate realizes concrete tracer classes, resolves parameter defaults and overrides, discovers required drivers, builds process-local topology, and compiles one static tracer equation for every realized tracer.

```julia
using Agate.Configuration: Population, Pool
using Agate.Processes: ModelDefinition

components = (
    N=Pool(:nitrogen),
    P=Population(; currency=:nitrogen, size_structure=[1.0]),
)

# Add named processes and their parameter definitions, then construct directly:
definition = ModelDefinition(; components, processes, parameters)
bgc = Agate.Construction.construct(definition)
```

The runnable examples develop this pattern in stages:

- [Franks et al. (1986) NPZ structure](@ref franks_npz) introduces pools, populations, growth, grazing, mortality, semantic parameters, and direct construction.
- [Detritus and bacterioplankton](@ref detritus_bacteria) shows structured material pools, heterotrophic POM consumption, living bacterial prey, and a reusable temperature factor.
- [Mixotrophy from process participation](@ref mixotrophy) shows one ordinary population participating in both growth and grazing.

For an advanced multi-currency example, `AgateRegistry.jl/examples/darwin_multicurrency.jl` inspects the DARWIN definition with Geider light limitation, Liebig DIN/PO4 limitation, fixed C-reference stoichiometry, and DOM/POM routing.

## Components describe structure

A component does not permanently declare itself to be a producer, consumer, prey item, mixotroph, or bacterioplankton. Those roles emerge from named process participation.

```julia
M = Population(; currency=:nitrogen, size_structure=[4.0])
```

The same `M` can appear as the population in a growth process and as the consumer in a grazing process. Likewise, POM remains a material [`Pool`](@ref); bacterioplankton that consume it remain living [`Population`](@ref) state and can separately participate as grazing prey.

## Named factors compose process rates

Growth and consumption can carry named multiplicative factors. Factor names are scientific identities, so declaration order does not change the model.

```julia
growth = Growth(;
    population=:P,
    factors=(
        light=Light(:smith; driver=:PAR),
        nutrients=NutrientResponse(:monod; resource=:N),
    ),
)
```

A formulation can still contain its own internal mathematics. For example, `Nutrients(:liebig)` combines its resource responses by Liebig limitation while remaining one named factor in the parent process.

## Parameters bind to scientific requirements

Each formulation declares semantic requirements. [`ParameterProvision`](@ref) connects those requirements to model-level parameter names, while [`ParameterDefinition`](@ref) supplies constructor-time defaults.

```julia
ParameterDefinition(
    ParameterSpec(
        :maximum_growth_rate,
        :vector;
        axes=:plankton,
        provides=ParameterProvision(
            :growth_P,
            (:factors, :light),
            :smith,
            :maximum_rate,
        ),
    ),
    FillDefault(2 / 86400),
)
```

[`DerivedDefault`](@ref) expresses defaults calculated from other parameters, including palatability and assimilation matrices. Explicit parameter overrides remain ordinary construction inputs:

```julia
bgc = construct(
    definition;
    parameter_overrides=(maximum_growth_rate=[3e-5, 0.0],),
)
```

## Construction is setup-time compilation

The authoring graph is resolved before runtime:

```text
components + named processes + parameter definitions
                         |
                         v
                    normalization
                         |
                         v
          topology + semantic parameter bindings
                         |
                         v
              typed process contributions
                         |
                         v
              static setup-time lowering
                         |
                         v
          one compiled equation per tracer
```

Runtime tracer kernels therefore evaluate concrete compiled terms rather than interpreting a dynamic process graph.
