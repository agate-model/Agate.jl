# Examples

This directory contains runnable scripts that demonstrate how to build and use
Agate biogeochemistry models.

## Running an example

From the repository root:

```bash
julia --project=. examples/box_model_factories.jl
```

Most examples assume you have instantiated the repo environment first:

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Dependencies

Many examples are written to mirror the documentation pages and therefore depend
on the Oceananigans/OceanBioME stack for time integration and light forcing. In
addition, plotting examples use CairoMakie.

- `box_model_factories.jl`: 0D box model walkthrough for factory construction
  and parameter overrides.
- `box_model_add_heterotroph.jl`: demonstrates extending a factory model with a
  new group and new equations.
- `1D_column.jl`: 1D water-column example (Oceananigans + OceanBioME).
