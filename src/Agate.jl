module Agate

include("Library/Library.jl")
include("Models/Models.jl")
include("Constructors/Constructors.jl")

using .Library
using .Models
using .Constructors

export compute_allometric_parameters
export define_tracer_functions
export construct_size_structured_NPZD

end # module
