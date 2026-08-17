# # [Exporting a model recipe] (@id export_model_recipe_example)

# This example shows how to save an authored Agate model definition and replay it later.
# Recipes preserve the scientific construction inputs while keeping the runtime grid explicit.

# ## Loading dependencies

using Agate.Construction: export_recipe, import_recipe
using Agate.Models: NiPiZD
using OceanBioME: BoxModelGrid

nothing #hide

# ## Construct a model and recipe

grid = BoxModelGrid()

bgc, recipe = NiPiZD.construct_with_recipe(;
    grid,
    size_structure=(;
        phytoplankton=(P=(n=3, min_esd=1.0, max_esd=10.0, splitting=:log_splitting),),
        zooplankton=(Z=(n=2, min_esd=20.0, max_esd=100.0, splitting=:linear_splitting),),
    ),
)

nothing #hide

# ## Export the recipe

recipe_path = tempname() * ".json"
export_recipe(recipe_path, recipe)

println("Wrote model recipe to ", recipe_path)

nothing #hide

# ## Reload and replay

loaded = import_recipe(recipe_path)
replayed = NiPiZD.construct(loaded; grid)

println(typeof(bgc))
println(typeof(replayed))

nothing #hide
