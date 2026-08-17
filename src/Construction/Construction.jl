"""Model construction utilities for factory-defined biogeochemistry models."""
module Construction

import Oceananigans
using ..Configuration: PFTSpecification

export construct_factory, construct_factory_with_manifest
export capture_model_recipe, replay_factory, resolve_construction_scalar_type
export recipe_family, recipe_factory
export ModelRecipe, ModelManifest
export encode_recipe, decode_recipe, export_recipe, import_recipe
export PFTSpecification

include("recipe.jl")
include("recipe_serialization.jl")
include("generator.jl")
include("construct.jl")

end
