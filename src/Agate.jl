module Agate

include("Factories/Factories.jl")
include("Equations/Equations.jl")
include("Library/Library.jl")
include("Configuration/Configuration.jl")
include("Runtime/Runtime.jl")
include("Construction/Construction.jl")
include("Models/Models.jl")
include("Introspection.jl")

using .Library
using .Factories
using .Configuration
using .Runtime
using .Models
using .Equations
using .Construction

export Library
export Models
export Factories
export Configuration
export Runtime
export Equations

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
