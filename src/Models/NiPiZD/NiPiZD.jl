module NiPiZD

include("Tracers.jl")
include("Constructor.jl")

using .Tracers
using .Constructor

export construct, instantiate

export DEFAULT_PHYTO_ARGS,
    DEFAULT_PHYTO_GEIDER_ARGS, DEFAULT_ZOO_ARGS, DEFAULT_INTERACTION_ARGS, DEFAULT_BGC_ARGS

end # module
