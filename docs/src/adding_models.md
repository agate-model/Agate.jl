# Adding a model

This page is for developers who want to integrate a new ecosystem model into Agate.
The intent is to keep the required surface area small and to make GPU support via parametric floating-point types (`FT`) straightforward.

## Checklist

A new model typically needs:

 1. A model module under `src/Models/<ModelName>/`.
 2. A small public constructor `Agate.Models.<ModelName>.construct` that forwards to `Agate.Constructor.construct_factory(factory; ...)`.
 3. A `parameter_directory(factory)` describing required parameter names and shapes.
 4. `default_parameters(factory, ctx, FT)` that returns defaults using the community context.
 5. A set of compiled dynamics / tracers that consume the parameters and update tendencies.
 6. Tests that exercise defaults, overrides, and GPU smoke (when available).

## Recommended structure

### 1) Define a factory

Define a concrete factory type for your model (for example `struct MyModelFactory <: AbstractBGCFactory end`) and implement the factory hooks needed by the generic constructor.

### 2) Declare parameters with a directory

Implement `parameter_directory(factory)` returning a `NamedTuple` of `ParameterSpec`s. Use:

  - `shape = :scalar | :vector | :matrix`
  - `axes = (:consumer, :prey)` when a matrix is consumer-by-prey

This enables early validation and provides the metadata used for interaction normalization.

### 3) Define consumer/prey roles

If your model uses role-aware consumer-by-prey interactions (e.g. predators only in a subset of groups), implement `default_roles(factory)`:

  - `default_roles(factory) = (consumers = (:Z, ...), prey = (:P, :D, ...))`

Return `nothing` for either field to include all classes on that axis. Overlap is allowed.

If you do not implement `default_roles`, Agate assumes all groups are both consumers and prey.

### 4) Provide defaults

Implement `default_parameters(factory, ctx, FT)` and return a `NamedTuple` keyed by the directory keys.

For axis-tagged interaction matrices, return rectangular defaults of size `(length(ctx.consumer_indices), length(ctx.prey_indices))`.

### 5) Dynamics should consume rectangular interactions

At runtime, the canonical matrices live in:

  - `p.interactions.palatability`
  - `p.interactions.assimilation`

Use the axis maps from `p.interactions` to map between consumer/prey indices and global plankton indices.

If your dynamics follow the consumer-by-prey pattern, prefer calling helpers from `Agate.Models.PredationSums` rather than re-implementing index plumbing.

### 6) Tests

Add at least:

  - constructor defaults
  - parameter overrides (scalar, vector, matrix)
  - interaction overrides (rectangular and provider-function forms)
  - a small GPU smoke test (gated) that checks rectangular matrices move to the device via `p.interactions.*`
