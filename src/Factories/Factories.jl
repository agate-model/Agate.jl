"""Registered model-family interfaces and parameter metadata used during construction."""
module Factories

export AbstractBGCFactory

"""Registration token for named biogeochemical model families."""
abstract type AbstractBGCFactory end

include("parameter_directory.jl")
include("interface.jl")

end
