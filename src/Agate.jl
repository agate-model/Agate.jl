"""Composable process-based plankton biogeochemistry for Oceananigans."""
module Agate

include("compat.jl")
include("ModelFamilies/ModelFamilies.jl")
include("Library/Library.jl")
include("Components/Components.jl")
include("Parameters/Parameters.jl")
include("Processes/Processes.jl")
include("Runtime/Runtime.jl")
include("Compilation/Compilation.jl")
include("Diagnostics/Diagnostics.jl")
include("Construction/Construction.jl")
include("Models/Models.jl")
include("Introspection.jl")

using .Processes: ModelDefinition

export Library
export Models
export ModelFamilies
export Parameters
export Components
export Processes
export Runtime
export Diagnostics
export Construction
export Introspection
export ModelDefinition

end # module
