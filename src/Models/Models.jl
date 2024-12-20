module Models

include("Biogeochemistry.jl")
include("Parameters.jl")

using .Biogeochemistry
using .Parameters

export compute_allometric_parameters
export define_tracer_functions

end # module
