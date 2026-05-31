# # [Exporting a model setup] (@id export_model_setup_example)

# This example shows how to save the inputs used to construct an Agate model.
# Saving the setup makes it easier to recreate the same model later, or to share
# the model configuration with someone else.

# ## Loading dependencies
# The example uses Agate.jl to construct the model and OceanBioME.jl to provide a
# simple box-model grid. Agate.Manifests provides the functions used to export and reload 
# the model setup.

using Agate
using Agate.Manifests: construct_from_manifest, export_manifest
using Agate.Models: NiPiZD
using OceanBioME: BoxModelGrid

nothing #hide

# ## Construct a model and setup record

# `NiPiZD.construct_with_manifest` accepts the same keywords as `NiPiZD.construct`.
# In addition to the biogeochemistry object, it returns a setup record describing
# the model recipe and constructor inputs. Agate stores this setup as a manifest,
# which is a JSON-compatible dictionary.

grid = BoxModelGrid()

bgc, manifest = NiPiZD.construct_with_manifest(;
    grid,
    phyto_size_structure=(n=3, min_esd=1.0, max_esd=10.0, splitting=:log_splitting),
    zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:linear_splitting),
)

nothing #hide

# ## Export the manifest

# `export_manifest` writes the setup to a JSON file. The exported file records the
# biological model configuration, but not the physical grid, so the grid remains
# explicit when the setup is loaded again.

manifest_path = tempname() * ".json"
export_manifest(manifest_path, manifest)

println("Wrote model manifest to ", manifest_path)

nothing #hide

# ## Reload the manifest

# `construct_from_manifest` reads the exported manifest and reconstructs the Agate
# biogeochemistry object. The grid is supplied separately because it describes the
# physical domain rather than the biological model setup.

reloaded_bgc = construct_from_manifest(manifest_path; grid)

println(typeof(bgc))
println(typeof(reloaded_bgc))

nothing #hide
