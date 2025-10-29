module Models

include("Parameters.jl")
include("DARWIN/DARWIN.jl")
include("NiPiZD/NiPiZD.jl")
include("NiPiZiHD/NiPiZiHD.jl")

using .Parameters
using .DARWIN
using .NiPiZD
using .NiPiZiHD

export compute_allometric_parameters, create_size_structured_params
export define_tracer_functions

end # module
