module Agate

include("box_model.jl")
include("functions.jl")
include("Library/Library.jl")
include("Models/Models.jl")
using .Models


export run_npzd_boxmodel
export create_bgc_struct, add_bgc_methods


end # module
