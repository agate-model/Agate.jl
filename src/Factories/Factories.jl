"""Factory interfaces.

This module contains the shared base types and hook points used by model families to
expose construction-time defaults.

Factories are intentionally explicit: they declare defaults for community structure and
compiled-equation builders, and optionally provide group-level hooks.

The parameter directory (`parameter_definitions`, `ParameterSpec`, etc.) also lives here
because it is factory-facing metadata used during construction-time validation.
"""
module Factories

export AbstractBGCFactory

"""Abstract supertype for biogeochemical model factories."""
abstract type AbstractBGCFactory end

# Factory-facing parameter metadata + defaults.
include("parameter_directory.jl")

# Default builders + optional hook points.
include("interface.jl")

end # module Factories
