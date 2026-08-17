"""Model construction utilities for factory-defined biogeochemistry models."""
module Construction

import Oceananigans
using ..Configuration: PFTSpecification

export construct_factory
export construct_factory_with_realization
export ModelRecipe, ModelRealization
export PFTSpecification

include("recipe.jl")
include("generator.jl")
include("construct.jl")

end
