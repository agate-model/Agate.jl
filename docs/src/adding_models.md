# Adding a model

This page is for developers who want to integrate a new ecosystem model into Agate.

## Checklist

A new model typically needs:

 1. A model module under `src/Models/<ModelName>/`.
 2. A small public constructor `Agate.Models.<ModelName>.construct` that forwards to `Agate.Construction.construct_factory(factory; ...)`.
 3. A single-source `parameter_definitions(factory)` that pairs each `ParameterSpec` with a default provider.
 4. Optionally, `derived_matrix_specs(factory)` (and `derivation_deps` for the derivation functions) for derived interaction matrices.
 5. A set of compiled dynamics / tracers that consume the parameters and update tendencies.
 6. Tests that exercise defaults and overrides.

## Recommended structure

### 1) Define a factory

Define a concrete factory type for your model (for example `struct MyModelFactory <: AbstractBGCFactory end`) and implement the factory hooks needed by the generic constructor.

At minimum you should implement:

  - `default_plankton_dynamics(factory)`
  - `default_biogeochem_dynamics(factory)`
  - `default_community(factory)`
  - `parameter_definitions(factory)` (spec + default in one place)
  - `parameter_directory(factory)` (derived automatically; can be overloaded if needed)

Canonical group ordering is inferred from the *explicit* ordering of the `community::NamedTuple` passed to the generic constructor.

### 2) Declare parameters with `parameter_definitions`

Implement `parameter_definitions(factory)` returning a tuple of `ParameterDefinition` entries (a `ParameterSpec` paired with a `DefaultProvider`).
Example:

```julia
using Agate.Utils: ParameterSpec, ParameterDefinition, ConstDefault, FillDefault, NoDefault

parameter_definitions(::MyModelFactory) = (
    ParameterDefinition(
        ParameterSpec(:remineralization, :scalar; doc="Detritus remineralization rate."),
        ConstDefault(0.18 / day),
    ),
    ParameterDefinition(
        ParameterSpec(:maximum_growth_rate, :vector; doc="Per-class growth rate."),
        FillDefault(2.0 / day),
    ),
    ParameterDefinition(
        ParameterSpec(:palatability_matrix, :matrix; axes=(:consumer, :prey), doc="Consumer preference."),
        NoDefault(),  # derived later
    ),
)
```
 Use:

  - `shape = :scalar | :vector | :matrix`
  - `axes = (:consumer, :prey)` when a matrix is consumer-by-prey

This enables early validation and provides the metadata used for interaction normalization.

### 3) Define consumer/prey roles

If your model uses role-aware consumer-by-prey interactions (e.g. predators only in a subset of groups), define a `roles` `NamedTuple` and pass it as `interaction_roles` to `Agate.Construction.construct_factory` (or set it as a default in your public constructor / model variant):

  - `roles = (consumers = (:Z, ...), prey = (:P, :D, ...))`

Return `nothing` for either field to include all classes on that axis. Overlap is allowed.

If you do not pass `interaction_roles`, Agate assumes all classes are both consumers and prey.

### 4) Default providers and derived matrices

Defaults are declared alongside parameter metadata in `parameter_definitions(factory)` by choosing an explicit `DefaultProvider`:

- `ConstDefault(x)` for scalar literals,
- `FillDefault(x)` for uniform vectors or matrices,
- `DiameterIndexedVectorDefault(x, :indices_field; default=...)` for vectors defined over a subset of diameter-indexed classes,
- `NoDefault()` when a parameter has no direct default (for example matrices derived later from trait vectors).

If your model derives interaction matrices from trait vectors, register the derivations with:

```julia
using Agate.Utils: derived_matrix_specs, derivation_deps

derived_matrix_specs(::MyModelFactory) = (;
    palatability_matrix = derive_palatability_matrix,
)

derivation_deps(::typeof(derive_palatability_matrix)) =
    (:optimum_predator_prey_ratio, :specificity, :protection)
```

The constructor recomputes a derived matrix when it is not explicitly overridden and at least one dependency key is explicitly overridden.


### 5) Dynamics should consume rectangular interactions

At runtime, the canonical matrices live in:

  - `p.interactions.palatability`
  - `p.interactions.assimilation`

Use the axis maps from `p.interactions` to map between consumer/prey indices and global plankton indices.

If your dynamics follow the consumer-by-prey pattern, prefer calling helpers from `Agate.Models.PredationSums` rather than re-implementing index plumbing.

### Optional: group-level hooks (OceanBioME-style)

Agate supports *group-level* dispatch via `Val{:Group}` for model-specific behaviour that should be chosen at compile time (GPU friendly).

Define methods in `Agate.Factories` such as:

```julia
using Agate.Factories: sinking_velocity, grazing_kernel

# Sinking speed (positive downward) for plankton group G and within-group class id
@inline sinking_velocity(model_or_factory, ::Val{:P}, class_id::Int, ctx, params) = nothing

# Grazing kernel selector for a predator/prey group pair (optional)
# Implementations should return a callable `k(prey, predator)`
@inline grazing_kernel(
    model_or_factory,
    ::Val{:Z},
    ::Val{:P},
    predator_class::Int,
    prey_class::Int,
    ctx,
    params,
) = my_kernel
```

These hooks are optional; if you do not implement them, Agate will use the default behaviour for your model.

### 6) Tests

Add at least:

  - constructor defaults
  - parameter overrides (scalar, vector, matrix)
  - interaction overrides (rectangular and provider-function forms)
  - a small GPU smoke test (gated) that checks rectangular matrices move to the device via `p.interactions.*`
