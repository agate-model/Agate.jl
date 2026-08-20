"""Model construction utilities for authored and registered biogeochemistry models."""
module Construction

import Oceananigans
using ..Configuration: PFTSpecification

export construct, construct_factory, construct_factory_plus_manifest
export ModelRecipe, ProcessModelRecipe, ModelManifest
export capture_process_model_recipe
export encode_recipe, decode_recipe, export_recipe, import_recipe
export PFTSpecification

include("recipe.jl")
include("recipe_serialization.jl")
include("recipe_provenance.jl")
include("generator.jl")
include("construct.jl")

end
