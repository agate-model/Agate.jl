# # [Exporting a model setup] (@id export_model_setup_example)

# This example shows how to construct a NiPiZD model, export the setup used to
# build it, and reload that setup later to reconstruct the model.

using Agate
using Agate.Manifests: construct_from_manifest, export_manifest
using Agate.Models: NiPiZD
using OceanBioME: BoxModelGrid

nothing #hide

# ## Construct with a setup record

# `NiPiZD.construct_with_manifest` accepts the same keywords as `NiPiZD.construct`,
# but also returns a JSON-compatible setup record. Agate represents this setup as
# a manifest.

grid = BoxModelGrid()

bgc, setup = NiPiZD.construct_with_manifest(;
    grid=grid,
    phyto_size_structure=(n=3, min_esd=1.0, max_esd=10.0, splitting=:log_splitting),
    zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:linear_splitting),
)

println(setup["model"]["id"])
println(setup["recipe"]["family"])
println(join(setup["resolved"]["tracers"], ", "))

nothing #hide

# ## Export

# Exporting writes the setup record to a JSON file. Construction does not write
# files automatically.

setup_path = tempname() * ".json"
export_manifest(setup_path, setup)

println("Wrote model setup to ", setup_path)

nothing #hide

# ## Import

# The exported setup can be loaded later to reconstruct the Agate biogeochemistry
# object. The physical grid is not encoded in the setup, so it is supplied when
# replaying it.

reloaded_bgc = construct_from_manifest(setup_path; grid)

println(typeof(bgc))
println(typeof(reloaded_bgc))

nothing #hide
