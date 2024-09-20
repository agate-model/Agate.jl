module Agate

include("functions.jl")
include("Growth/Growth.jl")
include("Models/Models.jl")
using .Models

export DynamicBGC, create_bgc_model

end # module
