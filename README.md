# Agate.jl

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/agate-model/Agate.jl/blob/main/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://agate-model.github.io/Agate.jl/dev/)
[![Build Status](https://github.com/agate-model/Agate.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/agate-model/Agate.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Aquatic GCM-Agnostic Tunable Ecosystems

**Agate.jl builds flexible and composable aquatic ecosystem models for [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) and [OceanBioME.jl](https://github.com/OceanBioME/OceanBioME.jl).**

## Installation

Agate is registered in the Julia General registry and requires Julia 1.10 or later.

Install it from the Julia package manager:

```julia
julia> ]
pkg> add Agate
```

or:

```julia
using Pkg
Pkg.add("Agate")
```

Then load Agate with:

```julia
using Agate
```

## Examples

The documentation includes examples covering:

* [Size-structured plankton communities](https://agate-model.github.io/Agate.jl/dev/generated/size_structure/)
* [Plankton functional types with different light strategies](https://agate-model.github.io/Agate.jl/dev/generated/named_pfts/)
* [Allometric parameter relationships](https://agate-model.github.io/Agate.jl/dev/generated/allometric_relationships/)
* [Predator-prey palatability](https://agate-model.github.io/Agate.jl/dev/generated/predator_prey_palatability/)
* [Models with mixotrophy](https://agate-model.github.io/Agate.jl/dev/generated/mixotrophy/)
* [Models with bacterioplankton and detritus consumption](https://agate-model.github.io/Agate.jl/dev/generated/detritus_bacteria/)
* [Forward-mode automatic differentiation](https://agate-model.github.io/Agate.jl/dev/generated/forward_mode_ad_nipizd_sensitivity/)
* [Reverse-mode automatic differentiation](https://agate-model.github.io/Agate.jl/dev/generated/reverse_mode_ad_nipizd_sensitivity/)
* [Exporting and replaying model definitions](https://agate-model.github.io/Agate.jl/dev/generated/export_model_recipe/)

## Documentation

The full documentation is available at:

**[agate-model.github.io/Agate.jl/dev/](https://agate-model.github.io/Agate.jl/dev/)**

Useful starting points include:

* [Quick start](https://agate-model.github.io/Agate.jl/dev/quick_start/)
* [Agate.jl-NiPiZD](https://agate-model.github.io/Agate.jl/dev/nipizd/)
* [Defining a model with mixotrophy](https://agate-model.github.io/Agate.jl/dev/generated/mixotrophy/)
* [API reference](https://agate-model.github.io/Agate.jl/dev/api/)

## Development and contributing

Agate is under active development with breaking changes still expected.

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup, testing, formatting, and contribution guidelines.
