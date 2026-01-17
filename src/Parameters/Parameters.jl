"""Agate.Parameters

Parameter registries and runtime resolution.

The parameter registry (`ParamRegistry`) is **CPU-side metadata** that declares:

- parameter names
- shapes (`:scalar`, `:vector`, `:matrix`)
- value kind (`:real` or `:bool`)
- documentation strings
- optional default providers

Runtime resolution happens during `Agate.Constructor.construct` and produces a
runtime parameter bundle containing only numeric scalars and dense arrays.
The runtime bundle is `Adapt.jl`-compatible so Oceananigans can adapt it to GPU.

Design notes
------------
- The registry is *explicit* and shape-driven: resolution never infers shapes from values.
- Vector parameters are strict-by-default. Group-level vectors use `GroupVec`.
- A small number of **built-in** derived matrices (currently palatability and assimilation)
  are computed during resolution when not explicitly provided.
"""
module Parameters

using Logging

using ..Library.Allometry: AbstractParamDef, resolve_param

using ..Utils.Specifications: ModelSpecification

import Base: show

export ParamSpec, ParamRegistry
export GroupVec
export scalar_param, vector_param, matrix_param
export parameter_registry, parameter_directory
export resolve_runtime_parameters
export update_registry

include("types.jl")
include("normalize.jl")
include("registry_ops.jl")
include("pretty.jl")
include("resolve.jl")

end # module Parameters
