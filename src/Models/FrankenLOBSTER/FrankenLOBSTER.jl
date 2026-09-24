"""Canonical FrankenLOBSTER model family and OceanBioME plankton integration boundary."""
module FrankenLOBSTER

include("nutrients.jl")
include("definition.jl")
include("parameters.jl")
include("interface.jl")
include("construction.jl")

export construct, construct_plus_recipe

end # module
