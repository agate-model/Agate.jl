module Configuration

using ..Parameters: ParameterSpec, parameter_directory, parameter_spec
import Adapt

export Population, Pool
export currency, states, state_currency, size_structure
export ComponentLayout, realize_components, realize_component_groups
export component_classes, component_state_tracers, component_tracers
export component_indices, state_tracers, state_indices
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
