module Agate

include("Library/Library.jl")
include("Models/Models.jl")

using .Library
using .Models

export compute_allometric_parameters
export define_tracer_functions

end # module
