module Agate

include("Library/Library.jl")
include("Models/Models.jl")

include("Models/Parameters.jl")

using .Models

export add_plankton_tracers
export define_tracer_functions

end # module
