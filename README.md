# Agate.jl &middot; [![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/agate-model/Agate.jl/blob/main/LICENSE) [![Documentation](https://img.shields.io/badge/docs-dev-blue)](https://agate-model.github.io/Agate.jl/dev/) [![Build Status](https://github.com/nanophyto/AGATE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/agate-model/Agate.jl/actions/workflows/CI.yml?query=branch%3Amain)
## Aquatic Gcm-Agnostic Tunable Ecosystems 

A Julia library to build flexible and composable aquatic ecosystems.


## Documentation

- [**DEV**](https://agate-model.github.io/Agate.jl/dev/) &mdash; *documentation of the in-development version.*

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

Which is (I think) equivalent to running the below in the Julia REPL (`]` enters Julia package manager mode):

```julia
]
activate Agate
instantiate
```

You can then use the package interactively, in the terminal:

```bash
julia --project=.
using Agate
```

or in a notebook:

```julia
using Pkg
Pkg.activate("Agate")

using Agate
```
