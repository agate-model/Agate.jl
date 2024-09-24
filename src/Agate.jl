module Agate

include("functions.jl")
include("Growth/Growth.jl")
include("Models/Models.jl")
using .Models

export create_struct, construct_bgc_model


end # module
