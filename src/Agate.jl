module Agate

include("Equations/Equations.jl")
include("Library/Library.jl")
include("Utils/Utils.jl")
include("Parameters/Parameters.jl")
include("Models/Models.jl")
include("Constructor/Constructor.jl")

using .Library
using .Utils
using .Parameters
using .Equations
using .Models
using .Constructor

export Library
export Models
export Utils
export Constructor
export Parameters
export Equations

# Re-export primary user-facing API
export construct
export default_community
export parameter_registry, parameter_directory, update_registry, extend_registry
export update_community, extend_community, update_dynamics, extend_dynamics
export NiPiZDFactory, DarwinFactory

export define_tracer_functions

end # module
