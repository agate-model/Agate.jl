module Models

include("Biogeochemistry.jl")
include("Parameters.jl")
include("Tracers.jl")

using .Biogeochemistry
using .Parameters
using .Tracers

export compute_allometric_parameters, create_params_dict
export define_tracer_functions

end # module
