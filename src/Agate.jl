module Agate

include("box_model.jl")
include("Library/Library.jl")
include("Models/Models.jl")
using .Models

export create_box_model, run_box_model
export create_bgc_struct, add_bgc_methods

end # module
