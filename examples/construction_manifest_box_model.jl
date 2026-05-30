# # [Construction manifest box model] (@id construction_manifest_box_model_example)

# !!! info
#     This example uses [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/) and [OceanBioME.jl](https://oceanbiome.github.io/OceanBioME.jl/stable/).
#     We recommend familiarizing yourself with their user interface if you intend to change the physical model setup.

# In this example we construct a DARWIN model for a zero-dimensional box model and export a
# construction manifest. The manifest is kept out of the runtime object, so the `bgc` passed
# to OceanBioME.jl stays focused on simulation.

# ## Loading dependencies
# The example uses Agate.jl for ecosystem construction and manifest export.
# Oceananigans.jl and OceanBioME.jl provide the default box-model setup.

using Agate
using Agate.Library.Light
using Agate.Models: export_manifest
using Agate.Models: DARWIN
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans

nothing #hide

# ## Ecosystem model and construction manifest

# `DARWIN.construct_with_manifest` accepts the same keywords as `DARWIN.construct`, but returns
# two values:
# - `bgc`: the runtime biogeochemistry object used by OceanBioME.jl
# - `manifest`: a JSON-compatible reconstruction recipe with resolved model metadata

bgc, manifest = DARWIN.construct_with_manifest(;
    grid=BoxModelGrid(),
    phyto_size_structure=(n=3, min_esd=1.5, max_esd=20.0, splitting=:log_splitting),
    zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:log_splitting),
)
nothing #hide

# The manifest includes a recipe section for reconstructing the Agate biogeochemistry
# object, plus a resolved section for inspecting the model that was built. The recipe
# stores explicit plankton diameters so future reconstruction does not need to rehydrate
# size-splitting options.

println(manifest["model"]["id"])
println(manifest["recipe"]["constructor"])
println(join(manifest["resolved"]["tracers"], ", "))
nothing #hide

# ## Physical model

# Next, we combine the Agate.jl ecosystem model with OceanBioME.jl.
# The `BoxModelGrid` represents a well-mixed zero-dimensional water column.

light_attenuation = FunctionFieldPAR(; grid=BoxModelGrid())
nothing #hide

bgc_model = Biogeochemistry(bgc; light_attenuation)
nothing #hide

full_model = BoxModel(; biogeochemistry=bgc_model)
nothing #hide

# ## Exporting the manifest

# Construction does not write files automatically.
# Export is explicit, and writes the already-created manifest to JSON.

manifest_path = tempname() * ".json"
export_manifest(manifest_path, manifest)

println("Wrote construction manifest to ", manifest_path)
nothing #hide
