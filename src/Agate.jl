module Agate

include("simulate.jl")

include("Library/Library.jl")
include("Models/Models.jl")

using .Models

export define_tracer_functions
export run_simulation

end # module
