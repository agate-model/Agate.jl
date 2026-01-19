# Beginner guide {#beginner-guide}

This section is for readers who are **new to Julia**, but have used Python / MATLAB / R.
The goal is to get you productive quickly.

## A short mental model

Agate helps you assemble aquatic ecosystem models from a few explicit ingredients:

- A **model constructor** (e.g. NiPiZD or DARWIN).
- A **community specification** (how many plankton size-classes, their diameter ranges, and which groups exist).
- **Parameters** (scalars, vectors, and matrices). You can override defaults with `NamedTuple`s.
- **Interaction matrices** (palatability and assimilation) that couple consumers to prey.
- An **architecture** choice (CPU or GPU).

A model constructor returns a biogeochemistry object you can plug into an Oceananigans or OceanBioME simulation. See the examples for end-to-end runs.

## Julia translation layer (coming from Python)

You do not need to "learn Julia" to use Agate, but a few correspondences help:

- **Keyword arguments**: in Julia you call `f(x=1, y=2)` (similar to Python). You may also see `f(; x=1, y=2)`.
  The semicolon `;` separates *positional* arguments from *keyword* arguments. If there are no positional arguments, the `;` is optional.
  See the Julia manual: [Keyword Arguments](https://docs.julialang.org/en/v1/manual/functions/#Keyword-Arguments).
- **NamedTuple**: `(; a=1, b=2)` is like an immutable dict with fixed keys; Agate uses `NamedTuple`s for configuration and overrides.
  See: [NamedTuple](https://docs.julialang.org/en/v1/base/base/#NamedTuple).

If you're curious later, Julia's "multiple dispatch" is why `construct` can exist for many model modules.
See: [Methods](https://docs.julialang.org/en/v1/manual/methods/).

## Getting started

Start here:

- Read [Quick start](@ref "Quick start") for a complete working setup.
- Pick a model: [NiPiZD](@ref NiPiZD) or [DARWIN](@ref DARWIN).
- If you are using a manuscript or project-specific configuration, see [Variants](@ref "Variants").
- Clone an example and modify one thing at a time.

## Common customisations

### Override a parameter

Most constructors accept a `parameters = (; ...)` keyword. For example:

```julia
params = (; maximum_growth_rate = 0.8)
model = Agate.Models.NiPiZD.construct(; parameters=params)
```

### Override interaction matrices

For palatability and assimilation matrices, the recommended override is a **rectangular consumer-by-prey** matrix of size
`(n_consumer, n_prey)`.

You can also provide a function `f(ctx)` that returns the matrix. Here, `ctx` is an **interaction context** describing the community
(size classes, groups, and which indices are consumers vs prey). In practice, you mainly use it to size your matrix:

```julia
pal_provider = ctx -> begin
    n_consumer = length(ctx.consumer_indices)
    n_prey = length(ctx.prey_indices)
    return ones(n_consumer, n_prey)
end

model = Agate.Models.NiPiZD.construct(; interactions=(; palatability_matrix=pal_provider))
```

Advanced forms (full square matrices and group-block matrices) exist for convenience, but you can ignore them when you're starting.
See the reference: [Parameters and interaction matrices](@ref "Parameters and interaction matrices").

## When you are ready for deeper changes

If you want to add a new model, define new parameters, or change the compiled dynamics, jump to the [Developer guide](@ref developer-guide).
