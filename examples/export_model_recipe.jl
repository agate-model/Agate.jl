# # [Exporting a model recipe] (@id export_model_recipe)

# Recipes record versioned family construction inputs for export and replay.

using Agate.Construction: export_recipe, import_recipe, encode_recipe
using Agate.Models: NiPiZD
using OceanBioME: BoxModelGrid
using Pkg

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
# family), a SHA-256 content fingerprint, and package provenance. Git repository and commit
# information are included when the implementation matches a clean checkout.

println("Recipe schema: ", encode_recipe(recipe)["schema"])
println("Saved ", basename(recipe_path))

nothing #hide

# ## Preserve the Julia environment
# The recipe identifies the scientific construction inputs; the Julia `Manifest.toml` pins
# the versions and source revisions of Agate, any external model provider, and dependencies.
# Archive the `Project.toml` and `Manifest.toml` of the environment used to construct the model.

project_path = Base.active_project()
manifest_path = joinpath(dirname(project_path), "Manifest.toml")
isfile(manifest_path) || Pkg.instantiate()
isfile(manifest_path) || error("The active Julia environment has no Manifest.toml")
cp(project_path, joinpath(release_dir, "Project.toml"))
cp(manifest_path, joinpath(release_dir, "Manifest.toml"))

println("Release files: ", join(sort(readdir(release_dir)), ", "))

nothing #hide

# For a portable release, use a **consumer project** with Agate and any model-provider
# package installed from registered versions or pinned Git revisions. Development environments
# containing `Pkg.develop(path=...)` entries (including this documentation build) refer to
# local source paths and must be replaced before distributing the environment.
# In a fresh session, install the pinned dependencies before replay:
#
# ```bash
# julia --project=/path/to/release -e 'using Pkg; Pkg.instantiate()'
# ```

# ## Reload and replay

# Runtime choices such as the grid are supplied again when the recipe is replayed.

loaded = import_recipe(recipe_path)
replayed = NiPiZD.construct(loaded; grid)

println("Recipe preserved: ", loaded == recipe)
println("Parameters preserved: ", replayed.parameters == bgc.parameters)

nothing #hide
