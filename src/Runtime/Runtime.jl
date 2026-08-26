# Runtime and diagnostic utilities.
module Runtime

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry, required_biogeochemical_auxiliary_fields, required_biogeochemical_tracers
using ..Configuration:
    ModelLayout, PopulationStateRef, component_class_count, component_indices,
    component_state_tracers, state_indices

export TracerValues
export tendency_inputs
export parameterized
export ode_problem
export active_parameters

include("tendency_context.jl")
include("active_parameters.jl")
include("ode_problem.jl")
include("class_refs.jl")
include("tracer_accessors.jl")

end # module
