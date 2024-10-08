# Agate.jl

[![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/agate-model/Agate.jl/blob/main/LICENSE)
[![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://agate-model.github.io/Agate.jl/dev/)
[![Build Status](https://github.com/agate-model/AGATE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/agate-model/Agate.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Aquatic Gcm-Agnostic Tunable Ecosystems

A Julia library to build flexible and composable aquatic ecosystems.

## Documentation

  - [**DEV**](https://agate-model.github.io/Agate.jl/dev/) — *documentation of the in-development version.*

## Getting started

Download Julia by following instructions at https://julialang.org/downloads/.

Clone this repository and change your current working directory to this project:

```bash
https://github.com/agate-model/Agate.jl.git
cd Agate.jl
```

To activate the project:

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Which is equivalent to running the below in the Julia REPL (`]` enters Julia package manager mode):

```Julia
]
activate Agate
instantiate
```

You can then use the package interactively, in the terminal:

```bash
julia --project=.
```

To use the package in a Jupyter notebook run:

```Julia
using Pkg
Pkg.activate("<path to Agate.jl repo>")
```

## Development

We follow the [Blue](https://github.com/JuliaDiff/BlueStyle) style guide for Julia. To automatically format all Julia files in the project you can use the JuliaFormatter. Once you have installed it (`add JuliaFormatter`) run:

```Julia
using JuliaFormatter

format(".")
```

To update project dependencies:

```Julia
] add <package>
```

To run tests:

```Julia
] test
```
