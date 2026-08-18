"""Model construction utilities for factory-defined biogeochemistry models."""
module Construction

import Oceananigans
using ..Configuration: PFTSpecification

export construct_factory, construct_factory_plus_manifest
export ModelRecipe, ModelManifest
export encode_recipe, decode_recipe, export_recipe, import_recipe
export PFTSpecification

include("recipe.jl")
include("recipe_serialization.jl")
include("recipe_provenance.jl")
include("generator.jl")
include("construct.jl")

end
