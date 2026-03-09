"""Utilities for size-dependent traits and interaction matrices."""
module Allometry

export AbstractParamDef, ConstantParam, AllometricParam
export PowerLaw
export PalatabilityPreyParameters, PalatabilityPredatorParameters
export allometric_scaling_power
export allometric_palatability_unimodal, allometric_palatability_unimodal_protection

include("parameter_defs.jl")
include("interactions.jl")

end # module Allometry
