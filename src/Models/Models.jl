module Models

include("Biogeochemistry.jl")
include("Constructors.jl")
include("Parameters.jl")
include("Tracers.jl")

using .Biogeochemistry
using .Constructors
using .Parameters
using .Tracers

export compute_allometric_parameters
export define_tracer_functions
export construct_NPZD_instance

end # module
