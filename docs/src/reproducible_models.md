# Reproducible model definitions

Agate separates a scientific model definition from its runtime environment.

- A [`ModelRecipe`](@ref Agate.Construction.ModelRecipe) captures authored scientific construction inputs before parameter and interaction materialization.
- A [`ModelManifest`](@ref Agate.Construction.ModelManifest) records the complete resolved deterministic scientific state produced by construction.
- Runtime objects such as the grid and architecture remain explicit environment inputs and are not stored in a recipe.

## Capture and replay

For the built-in NiPiZD family, `construct_with_manifest` returns the model, its recipe, and its resolved manifest:

```julia
using Agate.Models: NiPiZD

bgc, recipe, manifest = NiPiZD.construct_with_manifest()
replayed, replayed_manifest = NiPiZD.construct_with_manifest(recipe)

replayed_manifest == manifest
```

Recipes preserve semantic size specifications, parameter definitions and authored overrides, interaction definitions and overrides, ecological and interaction roles, sinking choices, open-bottom behavior, and scalar precision. Partial parameter overrides and parameter laws therefore remain authored inputs rather than being replaced by resolved vectors.

## Typed JSON recipes

The recipe document schema is `agate.model_recipe.v1`. Recipe serialization is typed so that supported semantic objects such as parameter laws, size specifications, matrices, symbols, scalar precision, and non-finite floating-point values can be reconstructed without relying on unrestricted reflection.

```julia
using Agate.Models: NiPiZD
using Agate.Construction: export_recipe, import_recipe

_, recipe = NiPiZD.construct_with_recipe()
export_recipe("nipizd.json", recipe)
loaded = import_recipe("nipizd.json")

loaded == recipe
```

The reader rejects unsupported schemas, unknown fields, unsupported model-family identifiers, and unsupported typed recipe values rather than silently ignoring them.

## External model families

External packages can extend Agate's model-family dispatch and provide their own constructors while reusing the generic recipe, manifest, parameter-law, interaction, and runtime machinery. The package defining a model family must be loaded before Agate can decode that family's recipe identifier.
