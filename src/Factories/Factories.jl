"""Factory interfaces and parameter metadata used during model construction."""
module Factories

export AbstractBGCFactory

"""Abstract supertype for biogeochemical model factories."""
abstract type AbstractBGCFactory end

include("parameter_directory.jl")
include("interface.jl")

end
