module Agate

include("box_model.jl")
include("bgc_model.jl")
include("Library/Library.jl")
include("Models/Models.jl")
using .Models

export define_tracer_functions
export create_bgc_model
export create_box_model, run_box_model


end # module
