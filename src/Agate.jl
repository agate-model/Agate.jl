module Agate

include("Utils/Utils.jl")
include("Equations/Equations.jl")
include("Library/Library.jl")
include("Parameters/Parameters.jl")
include("FactoryInterface.jl")
include("Constructor/Constructor.jl")
include("Models/Models.jl")
include("Introspection.jl")

using .Library
using .Utils
using .Parameters
using .FactoryInterface
using .Equations
using .Models
using .Constructor

export Library
export Models
export Utils
export Parameters
export Equations
export FactoryInterface

# Public model modules.
const NiPiZD = Models.NiPiZD
const DARWIN = Models.DARWIN

export NiPiZD
export DARWIN

# Newcomer UX helpers
export tracer_names
export auxiliary_field_names
export parameter_names
export required_parameters
export model_summary
export describe

end # module
