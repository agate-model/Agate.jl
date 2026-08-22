module Configuration

using ..Parameters: parameter_directory

export Population, Pool, PopulationStateRef, population_state
export currency, states, state_currency, size_structure
export ComponentLayout, realize_components, realize_component_groups
export component_classes, component_state_tracers, component_tracers
export component_indices, state_tracers, state_tracer, state_indices
export component_diameters, component_class_count

export PFTSpecification

export build_plankton_community

export PalatabilityAllometric, AssimilationBinary

include("specifications.jl")
include("community.jl")
include("components.jl")
include("interactions_matrices.jl")
include("interactions_derivations.jl")

end # module
