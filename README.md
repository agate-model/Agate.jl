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
Pkg.develop(path="/path/to/Agate.jl")
```

or add it directly from Git:

```julia
using Pkg
Pkg.add(url="https://github.com/agate-model/Agate.jl.git")
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
  - `docs/src/generated/`: precomputed outputs (figures and small datasets) used by the docs build.
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
