module Configuration

using ..Factories: AbstractBGCFactory, ParameterSpec, parameter_spec
import Adapt

export PFTSpecification, pft_get, pft_has

export build_plankton_community
export AbstractDiameterSpecification, DiameterListSpecification, DiameterRangeSpecification
export diameter_specification
export validate_community_inputs, parse_community, CommunityContext

export axis_indices
export normalize_interaction_overrides
export InteractionMatrices
export finalize_interaction_parameters

export derivation_deps, derived_matrix_specs, resolve_derived_matrices
export derive_palatability_matrix_allometric, derive_assimilation_matrix_binary

include("specifications.jl")
include("community.jl")
include("interactions_matrices.jl")
include("interactions_derivations.jl")

end # module
