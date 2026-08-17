# # [Exporting a model recipe] (@id export_model_recipe_example)

# This example shows how to save an authored Agate model definition and replay it later.
# Recipes preserve the scientific construction inputs while keeping the runtime grid explicit.

# ## Loading dependencies

using Agate.Construction: export_recipe, import_recipe
using Agate.Models: NiPiZD
using OceanBioME: BoxModelGrid

nothing #hide

# ## Construct a model, recipe, and manifest

# `ModelRecipe` stores the scientific inputs used to construct the model. `ModelManifest`
# records the resolved scientific state produced from those inputs.

grid = BoxModelGrid()

bgc, recipe, manifest = NiPiZD.construct_with_manifest(;
    grid,
    size_structure=(;
        phytoplankton=(P=(n=3, min_esd=1.0, max_esd=10.0, splitting=:log_splitting),),
        zooplankton=(Z=(n=2, min_esd=20.0, max_esd=100.0, splitting=:linear_splitting),),
    ),
)

nothing #hide

# ## Export the recipe

# Recipes are the durable model definition. The resolved manifest stays in memory and can
# be used to check that replay produces the same scientific state.

recipe_path = tempname() * ".json"
export_recipe(recipe_path, recipe)

println("Wrote model recipe to ", recipe_path)

nothing #hide

# ## Reload and replay

# The grid is supplied again when replaying because runtime environment choices are not
# stored in the recipe.

loaded = import_recipe(recipe_path)
replayed, replayed_manifest = NiPiZD.construct_with_manifest(loaded; grid)

println("Recipe preserved: ", loaded == recipe)
println("Manifest preserved: ", replayed_manifest == manifest)
println(typeof(replayed))

nothing #hide
