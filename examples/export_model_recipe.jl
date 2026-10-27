# # [Exporting a model recipe] (@id export_model_recipe)

# Recipes record versioned family construction inputs for export and replay.

using Agate.Construction: export_recipe, import_recipe, encode_recipe, export_model_release
using Agate.Models: NiPiZD
using OceanBioME: BoxModelGrid

nothing #hide

# ## Construct and save
# `construct_plus_recipe()` returns the model together with the authored scientific recipe.
# The recipe can be exported to JSON and later reloaded to reconstruct the model.

grid = BoxModelGrid()

bgc, recipe = NiPiZD.construct_plus_recipe(;
    grid,
    size_structure=(;
        phytoplankton=(P=(n=3, min_esd=1.0, max_esd=10.0, spacing=:log),),
        zooplankton=(Z=(n=2, min_esd=20.0, max_esd=100.0, spacing=:linear),),
    ),
)

release_dir = mktempdir()
recipe_path = joinpath(release_dir, "recipe.json")
export_recipe(recipe_path, recipe)

# The JSON records the family definition version, declared ancestry (or `null` for a root
# family), a SHA-256 content fingerprint, and package provenance.

println("Recipe schema: ", encode_recipe(recipe)["schema"])
println("Saved ", basename(recipe_path))

nothing #hide

# ## Export a reproducible model release
# `export_model_release` saves the recipe together with the active Julia project's
# `Project.toml` and the manifest selected by Julia. It first checks that the
# recipe and loaded package sources match the active environment and that package
# sources are portable (with source-tree hashes for Git dependencies).
#
# Run the following from a Julia project with Agate and your model package installed:
#
# ```julia
# export_model_release("my-model-release", recipe)
# ```
#
# This documentation example runs from a development checkout, so it demonstrates
# recipe export and replay; the checked release export is for a portable Julia project.

# ## Reload and replay

# Runtime choices such as the grid are supplied again when the recipe is replayed.

loaded = import_recipe(recipe_path)
replayed = NiPiZD.construct(loaded; grid)

println("Recipe preserved: ", loaded == recipe)
println("Parameters preserved: ", replayed.parameters == bgc.parameters)

nothing #hide
