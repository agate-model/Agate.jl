module Configuration

using ..Factories: AbstractBGCFactory, ParameterSpec, parameter_directory, parameter_spec
import Adapt

export PFTSpecification

export build_plankton_community

export AbstractMatrixDeriver
export MatrixDefinition
export matrix_definitions
export derive_matrix, derivation_deps
export matrix_deriver_identifier, matrix_deriver_from_identifier
export PalatabilityAllometric, AssimilationBinary

include("specifications.jl")
include("community.jl")
include("interactions_matrices.jl")
include("interactions_derivations.jl")

end # module
