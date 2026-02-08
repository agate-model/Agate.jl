module Agate

include("Factories/Factories.jl")
include("Utils/Utils.jl")
include("Equations/Equations.jl")
include("Library/Library.jl")
include("Interface.jl")
include("Construction/Construction.jl")
include("Models/Models.jl")
include("Introspection.jl")

using .Library
using .Factories
using .Utils
using .Interface
using .Models
using .Equations
using .Construction

export Library
export Models
export Factories
export Utils
export Equations
export Interface

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
