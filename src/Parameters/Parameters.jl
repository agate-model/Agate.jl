"""Agate.Parameters

Parameter registries and resolution utilities.

The parameter registry is **CPU-only metadata**. Parameter *resolution* happens during
`construct` and produces a runtime parameter bundle containing only numeric scalars and
arrays. The returned runtime bundle is `Adapt.jl`-compatible, so it can be adapted to GPU
arrays by Oceananigans.

Provider system
---------------
Parameters are declared with an explicit `shape` (`:scalar`, `:vector`, or `:matrix`).
User inputs are normalized into canonical provider values at **boundary time**
(`ParamSpec`, `update_registry`, `extend_registry`, and interactions application).
Runtime resolution never guesses a shape from the value type.

Derived providers may declare dependencies via `deps(provider)::Tuple{Vararg{Symbol}}` so the
resolver can compute prerequisite parameters without hardcoded special cases.
"""
module Parameters

using Logging

using Agate.Library.Allometry: AbstractParamDef, resolve_param
using Agate.Library.Allometry:
    allometric_palatability_unimodal_protection,
    palatability_matrix_allometric,
    assimilation_efficiency_matrix_binary

using Agate.Utils.Specifications: PFTSpecification, pft_has, pft_get, ModelSpecification

import Base: show

export ParamSpec, ParamRegistry
export MatrixFn
export scalar_param, vector_param, matrix_param
export parameter_registry, parameter_directory
export resolve_runtime_parameters
export update_registry, extend_registry

include("types.jl")
include("normalize.jl")
include("registry_ops.jl")
include("pretty.jl")
include("defaults.jl")
include("resolve.jl")

end # module Parameters
