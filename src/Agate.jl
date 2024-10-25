module Agate

include("box_model.jl")
include("Library/Library.jl")
include("Models/Models.jl")
using .Models

export create_box_model, run_box_model
export define_tracer_functions

end # module
