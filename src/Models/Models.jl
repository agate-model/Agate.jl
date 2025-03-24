module Models

include("Parameters.jl")
include("DARWIN/DARWIN.jl")
include("NiPiZD/NiPiZD.jl")

using .Parameters
using .DARWIN
using .NiPiZD

export compute_allometric_parameters, create_params_dict
export define_tracer_functions

end # module
