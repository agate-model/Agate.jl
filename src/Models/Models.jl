module Models

include("Dynamic.jl")
include("Parameters.jl")

using .Dynamic
using .Parameters

export compute_darwin_parameters
export define_tracer_functions

end # module
