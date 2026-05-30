# # [Construction manifest box model] (@id construction_manifest_box_model_example)

# !!! info
#     This example uses [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/stable/) and [OceanBioME.jl](https://oceanbiome.github.io/OceanBioME.jl/stable/).
#     We recommend familiarizing yourself with their user interface if you intend to change the physical model setup.

# In this example we construct a default Agate.jl variant for use in a zero-dimensional box model.
# The construction workflow returns the runtime biogeochemistry object and a separate construction manifest.
# The manifest is kept out of the runtime object, so the `bgc` passed to OceanBioME.jl stays focused on simulation.

# ## Loading dependencies
# The example uses Agate.jl for ecosystem construction and manifest export.
# Oceananigans.jl and OceanBioME.jl provide the default box-model setup.

using Agate
using Agate.Library.Light
using Agate.Models: ModelId, variant, construct_with_manifest, export_manifest
using Agate.Models.DARWIN
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans

nothing #hide

# ## Ecosystem model and construction manifest

# First, we look up a built-in DARWIN variant and construct it with manifest generation enabled.
# `construct_with_manifest` returns two values:
# - `bgc`: the runtime biogeochemistry object used by OceanBioME.jl
# - `manifest`: a JSON-compatible record of the resolved construction state

id = ModelId(:DARWIN, :citation2026, :A)
spec = variant(id)

bgc, manifest = construct_with_manifest(spec; grid=BoxModelGrid())
nothing #hide

# The manifest records the model identity, Agate version, Julia version, resolved tracers,
# auxiliary fields, plankton diameters, and resolved parameter values.

println(manifest["model"]["id"])
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
