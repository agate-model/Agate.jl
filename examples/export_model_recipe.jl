# # [Exporting a model recipe] (@id export_model_recipe_example)

# Recipes store the authored scientific model definition so it can be exported and replayed.

using Agate.Construction: export_recipe, import_recipe
using Agate.Models: NiPiZD
using OceanBioME: BoxModelGrid

nothing #hide

# ## Construct and save

grid = BoxModelGrid()

bgc, recipe = NiPiZD.construct_plus_recipe(;
    grid,
    size_structure=(;
        phytoplankton=(P=(n=3, min_esd=1.0, max_esd=10.0, splitting=:log_splitting),),
        zooplankton=(Z=(n=2, min_esd=20.0, max_esd=100.0, splitting=:linear_splitting),),
    ),
)

recipe_path = tempname() * ".json"
export_recipe(recipe_path, recipe)

# The JSON also records a SHA-256 recipe fingerprint and package provenance. Git
# repository and commit information are included when the implementation matches a checkout.

println("Wrote model recipe to ", recipe_path)

nothing #hide

# ## Reload and replay

# Runtime choices such as the grid are supplied again when the recipe is replayed.

loaded = import_recipe(recipe_path)
replayed = NiPiZD.construct_from_recipe(loaded; grid)

println("Recipe preserved: ", loaded == recipe)
println("Parameters preserved: ", replayed.parameters == bgc.parameters)

nothing #hide
