"""Model construction utilities for factory-defined biogeochemistry models."""
module Construction

import Oceananigans
using ..Configuration: PFTSpecification

export construct_factory
export ModelRecipe, ModelRealization
export PFTSpecification

include("recipe.jl")
include("generator.jl")
include("construct.jl")

end
