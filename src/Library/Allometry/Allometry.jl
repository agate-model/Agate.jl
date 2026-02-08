"""Utilities for size-dependent traits and interaction matrices."""
module Allometry

export AbstractParamDef, ConstantParam, AllometricParam
export PowerLaw, resolve_param, cast_paramdef
export resolve_diameter_vector, resolve_diameter_indexed_vector
export PalatabilityPreyParameters, PalatabilityPredatorParameters
export allometric_scaling_power
export allometric_palatability_unimodal, allometric_palatability_unimodal_protection
export palatability_matrix_allometric_axes, assimilation_efficiency_matrix_binary_axes

include("parameter_defs.jl")
include("interactions.jl")

end # module Allometry
