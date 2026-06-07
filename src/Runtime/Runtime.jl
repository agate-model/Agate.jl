# Runtime and diagnostic utilities.
module Runtime

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

export TendencyContext
export TracerValues
export tendency_inputs
export parameterized

include("tendency_context.jl")
include("active_parameters.jl")
include("class_refs.jl")
include("tracer_accessors.jl")

end # module
