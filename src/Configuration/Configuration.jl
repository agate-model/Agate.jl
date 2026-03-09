module Configuration

using ..Factories: AbstractBGCFactory, ParameterSpec, parameter_spec
import Adapt

export PFTSpecification

export build_plankton_community

export AbstractMatrixDeriver
export MatrixDefinition
export matrix_definitions

include("specifications.jl")
include("community.jl")
include("interactions_matrices.jl")
include("interactions_derivations.jl")

end # module
