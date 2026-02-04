# Agate.jl

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/agate-model/Agate.jl/blob/main/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://agate-model.github.io/Agate.jl/dev/)
[![Build Status](https://github.com/agate-model/Agate.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/agate-model/Agate.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Aquatic Gcm-Agnostic Tunable Ecosystems

A Julia library to build flexible and composable aquatic ecosystems.

## Documentation

  - [**DEV**](https://agate-model.github.io/Agate.jl/dev/) — *documentation of the in-development version.*

## Getting started

Download Julia by following instructions at https://julialang.org/downloads/.

Clone this repository and change your current working directory to this project:

```bash
git clone https://github.com/agate-model/Agate.jl.git
cd Agate.jl
```

Instantiate the project environment (this may take a while the first time):

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

You can then use the package interactively, in the terminal:

```bash
julia --project=.
```

To run one of the scripts in `examples/`:

```bash
julia --project=. examples/box_model_factories.jl
```

To use the package from another Julia project, either `develop` a local checkout:

```julia
using Pkg
Pkg.develop(; path="/path/to/Agate.jl")
```

or add it directly from Git:

```julia
using Pkg
Pkg.add(; url="https://github.com/agate-model/Agate.jl.git")
```

In a Jupyter notebook, activate a project that has Agate in its environment:

```julia
using Pkg
Pkg.activate("/path/to/project")
```

## Repository layout

  - `src/`: package source code.

  - `test/`: test suite.
  - `examples/`: runnable scripts (many depend on Oceananigans/OceanBioME; see `examples/README.md`).
  - `docs/`: documentation build.
    
      + `docs/src/generated/`: precomputed outputs (figures and small datasets) used by the docs build.
  - `paper/`: scripts and container recipes used for paper/plot reproduction (see `paper/README.md`).

## Development

We follow the [Blue](https://github.com/JuliaDiff/BlueStyle) style guide for Julia. To automatically format all Julia files in the project you can use JuliaFormatter:

```julia
using Pkg
Pkg.add("JuliaFormatter")
using JuliaFormatter
format(".")
```

To update project dependencies in the REPL, press `]` and run:

```text
add <package>
```

To run tests:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

## Interaction matrices

Agate models expose two predator-by-prey interaction matrices:

  - `palatability_matrix` — preference of each consumer for each prey
  - `assimilation_matrix` — assimilation efficiency of each consumer on each prey

These matrices are **role-aware** and are stored canonically as
`(n_consumer, n_prey)` rectangular matrices.

All public model constructors accept overrides via the two keywords
`palatability_matrix=` and `assimilation_matrix=`. Each may be:

  - a rectangular `(n_consumer, n_prey)` matrix
  - a full `(n_total, n_total)` matrix (embedded/sliced automatically)
  - a provider function `(ctx) -> matrix` evaluated once at construction time

For advanced workflows, the construction context provides explicit axes:

```julia
pal = (ctx) -> begin
    nc = length(ctx.consumer_indices)
    np = length(ctx.prey_indices)
    M = zeros(ctx.FT, nc, np)
    # fill M in consumer-by-prey order
    return M
end

bgc = NiPiZD.construct(; palatability_matrix=pal)
```

If you have a small *group-by-group* block matrix, wrap it in
`Agate.Utils.GroupBlockMatrix(B)` to force expansion across all size classes.

## Derived matrices (trait overrides)

NiPiZD and DARWIN expose a small set of **interaction traits** (vectors) that
are used to derive default interaction matrices. If you override one of these
traits (and do not explicitly override the corresponding matrix), Agate
recomputes the matrix during construction.

Example (tighten palatability specificity for consumers):

```julia
n_phyto = 4
n_zoo = 2
n_total = n_phyto + n_zoo

bgc = NiPiZD.construct(; n_phyto, n_zoo, parameters=(; specificity=fill(0.15f0, n_total),))
```
