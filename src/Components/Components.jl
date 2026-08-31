"""Scientific components, diameter specifications, and realized model topology."""
module Components

export Plankton, Pool
export element, states, reference_state, variable_states, state_element, size_structure
export component_entities, component_state_tracers
export component_tracers, state_tracers, state_tracer
export component_diameters

include("diameters.jl")
include("component_types.jl")
include("layout.jl")

end # module
