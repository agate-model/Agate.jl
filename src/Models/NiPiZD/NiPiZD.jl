"""Bundled NiPiZD model family."""
module NiPiZD

include("definition.jl")
include("parameters.jl")
include("construction.jl")

export construct, construct_plus_recipe, construct_from_recipe

end # module
