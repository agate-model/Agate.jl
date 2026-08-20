module Configuration

using ..Parameters: ParameterSpec, parameter_directory, parameter_spec
import Adapt

export Population, Pool
export currency, size_structure, sinking
export ComponentLayout, realize_components, realize_component_groups
export component_tracers, component_indices, component_diameters, component_class_count
export realize_component_sinking

export PFTSpecification

export build_plankton_community

export PalatabilityAllometric, AssimilationBinary

include("specifications.jl")
include("community.jl")
include("components.jl")
include("interactions_matrices.jl")
include("interactions_derivations.jl")

end # module
