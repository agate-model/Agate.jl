# API reference

## Introspection

```@docs
Agate.Introspection.model_summary
Agate.Introspection.describe
Agate.Introspection.tracer_names
Agate.Introspection.auxiliary_field_names
Agate.Introspection.parameter_names
Agate.Introspection.plankton_groups
Agate.Introspection.plankton_tracers
Agate.Introspection.plankton_diameters
Agate.Introspection.nonplankton_tracers
Agate.Introspection.tracer_groups
Agate.Introspection.interaction_matrix
```

### Interaction matrix introspection

```julia
using Agate
using Agate.Introspection

bgc = Agate.Models.NiPiZD.construct()
pal = interaction_matrix(bgc, :palatability)

pal.rows     # consumer labels
pal.columns  # prey labels
pal.matrix   # consumer-by-prey matrix
```

## Construction API

### Recipes and low-level construction records

Model users normally work through the model-family constructors (`construct`,
`construct_plus_recipe`, and `construct_from_recipe`). The manifest API below is a
lower-level resolved-state record used for diagnostics, model-family development, and
semantic replay tests across explicit execution environments.

`ModelRecipe` uses the `agate.model_recipe.v2` schema and records scientific construction
semantics only. Fresh model inputs canonicalize to a recipe, and loaded recipes enter the
same `construct_factory(recipe; ...)` realization path. Scalar precision is resolved from
the execution environment (or an explicit `scalar_type` override) and is recorded on
`ModelManifest`, not on the recipe.

```@docs
Agate.Construction.ModelRecipe
Agate.Construction.ModelManifest
Agate.Construction.construct_factory
Agate.Construction.construct_factory_plus_manifest
Agate.Construction.encode_recipe
Agate.Construction.decode_recipe
Agate.Construction.export_recipe
Agate.Construction.import_recipe
```

### External model-family API

Extension hooks are imported explicitly by packages that define model families.

```@docs
Agate.Construction.capture_model_recipe
Agate.Construction.replay_factory
Agate.Construction.resolve_construction_scalar_type
Agate.Construction.recipe_family
Agate.Construction.recipe_factory
Agate.Factories.AbstractBGCFactory
Agate.Factories.parameter_definitions
Agate.Factories.parameter_directory
Agate.Factories.parameter_spec
Agate.Factories.ParameterSpec
Agate.Factories.ParameterDefinition
Agate.Factories.DefaultProvider
Agate.Factories.ConstDefault
Agate.Factories.NoDefault
Agate.Factories.FillDefault
Agate.Factories.DiameterIndexedVectorDefault
Agate.Factories.DiameterIndexedMaterialization
Agate.Factories.default_components
Agate.Factories.default_community
Agate.Factories.default_plankton_dynamics
Agate.Factories.default_biogeochem_dynamics
Agate.Configuration.AbstractMatrixDeriver
Agate.Configuration.MatrixDefinition
Agate.Configuration.matrix_definitions
Agate.Configuration.derive_matrix
Agate.Configuration.derivation_deps
```

### Components and configuration

```@docs
Agate.Configuration.Population
Agate.Configuration.Pool
Agate.Configuration.currency
Agate.Configuration.size_structure
Agate.Configuration.sinking
Agate.Configuration.ComponentLayout
Agate.Configuration.realize_components
```

### Plankton communities and interactions

```@docs
Agate.Configuration.PFTSpecification
Agate.Configuration.build_plankton_community
Agate.Configuration.PalatabilityAllometric
Agate.Configuration.AssimilationBinary
```
