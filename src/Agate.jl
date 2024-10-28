module Agate

include("bgc_model.jl")
include("simulate.jl")

include("Library/Library.jl")
include("Models/Models.jl")

using .Models

export define_tracer_functions
export create_bgc_model
export run_simulation

end # module
