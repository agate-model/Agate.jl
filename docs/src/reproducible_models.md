# Reproducible model definitions

Agate separates authored scientific definitions, resolved scientific state, and runtime environment inputs.

- A [`ModelRecipe`](@ref Agate.Construction.ModelRecipe) captures authored scientific construction inputs before parameter and interaction materialization.
- A [`ModelManifest`](@ref Agate.Construction.ModelManifest) records the complete resolved deterministic scientific state produced by construction.
- Runtime objects such as the grid and architecture remain explicit environment inputs and are not stored in a recipe.

Recipes are the durable serialized model definitions. Model manifests are resolved in-memory construction records and currently have no serialized format.

## Capture and replay

For NiPiZD, `construct_with_manifest` returns the model, its recipe, and its resolved manifest:

```julia
using Agate.Models: NiPiZD

bgc, recipe, manifest = NiPiZD.construct_with_manifest()
replayed, replayed_manifest = NiPiZD.construct_with_manifest(recipe)

replayed_manifest == manifest
```

Recipes preserve semantic size specifications, parameter definitions and authored overrides, interaction definitions and overrides, ecological and interaction roles, sinking choices, open-bottom behavior, and scalar precision. Partial parameter overrides and parameter laws therefore remain authored inputs rather than being replaced by resolved vectors.

## Typed JSON recipes

The recipe document schema is `agate.model_recipe.v1`. Recipe serialization is typed so that supported semantic objects such as parameter laws, size specifications, matrices, symbols, scalar precision, and non-finite floating-point values can be reconstructed without unrestricted reflection.

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

External packages can define model families while reusing Agate's generic community, parameter, interaction, recipe, manifest, and runtime machinery. The package defining a family must be loaded before Agate can decode that family's recipe identifier.

The supported extension contract is namespaced by responsibility.

### Factory definition

A family defines an [`AbstractBGCFactory`](@ref Agate.Factories.AbstractBGCFactory) subtype and implements:

- `Agate.Factories.parameter_definitions(factory)` for parameter metadata and constructor-time defaults;
- `Agate.Factories.default_community(factory)` for structural group and size defaults;
- `Agate.Factories.default_plankton_dynamics(factory)` for plankton tendency builders;
- `Agate.Factories.default_biogeochem_dynamics(factory)` for non-plankton tendency builders;
- `Agate.Configuration.matrix_definitions(factory)` when interaction matrices are derived from other parameters;
- `Agate.Construction.recipe_family(factory)` and `Agate.Construction.recipe_factory(::Val{family})` for stable recipe identity and decoding.

A minimal family skeleton is:

```julia
using Agate.Configuration: PFTSpecification
using Agate.Factories: AbstractBGCFactory
import Agate.Construction: recipe_family, recipe_factory
import Agate.Factories:
    parameter_definitions,
    default_community,
    default_plankton_dynamics,
    default_biogeochem_dynamics

struct ExampleFactory <: AbstractBGCFactory end

recipe_family(::ExampleFactory) = :Example
recipe_factory(::Val{:Example}) = ExampleFactory()

parameter_definitions(::ExampleFactory) = ()

default_community(::ExampleFactory) = (
    P=(; diameters=[2.0], pft=PFTSpecification()),
)

# `example_phytoplankton` is a model-specific dynamics builder.
default_plankton_dynamics(::ExampleFactory) = (P=example_phytoplankton,)
default_biogeochem_dynamics(::ExampleFactory) = (;)
```

The model package remains responsible for its scientific constructor surface: grouping user inputs into recipe inputs versus runtime environment inputs, declaring ecological/interaction/parameter roles, and selecting its tendency builders.

### Parameter definitions

The stable parameter-definition surface is:

- `Agate.Factories.ParameterSpec` and `ParameterDefinition`;
- `ConstDefault`, `NoDefault`, `FillDefault`, and `DiameterIndexedVectorDefault`;
- `DiameterIndexedMaterialization` for declaring where parameter-law overrides apply;
- `parameter_definitions`, `parameter_directory`, and `parameter_spec`.

Default generation and parameter-law materialization are separate semantics. A parameter's default provider does not determine whether or where a parameter law may be materialized.

### Community and interactions

External families may use:

- `Agate.Configuration.PFTSpecification`;
- `build_plankton_community` and `normalize_diameters`;
- `AbstractMatrixDeriver`, `MatrixDefinition`, and `matrix_definitions`;
- the built-in `PalatabilityAllometric` and `AssimilationBinary` derivers.

A custom matrix deriver implements `derive_matrix`. If the deriver is stored in a serialized recipe it must also implement a stable `matrix_deriver_identifier` and matching `matrix_deriver_from_identifier(::Val{...})` method.

### Recipe capture and construction

Model packages may rely on the following construction helpers:

- `Agate.Construction.capture_model_recipe`;
- `replay_factory`;
- `resolve_construction_scalar_type`;
- `construct_factory` and `construct_factory_with_manifest`;
- `ModelRecipe` and `ModelManifest`;
- `encode_recipe`, `decode_recipe`, `export_recipe`, and `import_recipe`.

`replay_factory(recipe)` reuses the parameter and interaction definitions captured by the recipe, so replay does not consult a model family's current defaults for those definitions.
