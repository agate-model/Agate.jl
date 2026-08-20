# API reference

## Direct model construction

[`ModelDefinition`](@ref) is the scientific model-definition container. Direct construction realizes intrinsic component size structure, resolves process-owned parameter requirements and required drivers, and compiles runtime tracer equations during setup.

```@docs
Agate.Construction.construct
Agate.Processes.ModelDefinition
```

### Components

```@docs
Agate.Configuration.Population
Agate.Configuration.Pool
Agate.Configuration.currency
Agate.Configuration.size_structure
Agate.Configuration.sinking
Agate.Configuration.ComponentLayout
Agate.Configuration.realize_components
Agate.Configuration.realize_component_groups
```

### Processes and factors

```@docs
Agate.Processes.AbstractProcess
Agate.Processes.AbstractFormulation
Agate.Processes.AbstractFactor
Agate.Processes.Growth
Agate.Processes.Light
Agate.Processes.NutrientResponse
Agate.Processes.Nutrients
Agate.Processes.Temperature
Agate.Processes.factors
Agate.Processes.Grazing
Agate.Processes.Consumption
Agate.Processes.Mortality
Agate.Processes.ProductRouting
Agate.Processes.FixedStoichiometry
Agate.Processes.Remineralization
```

### Parameter definitions

```@docs
Agate.Factories.ParameterSpec
Agate.Factories.ParameterProvision
Agate.Factories.ParameterDefinition
Agate.Factories.DefaultProvider
Agate.Factories.ConstDefault
Agate.Factories.DerivedDefault
Agate.Factories.derive_default
Agate.Factories.NoDefault
Agate.Factories.FillDefault
Agate.Factories.DiameterIndexedVectorDefault
Agate.Factories.DiameterIndexedMaterialization
Agate.Configuration.PalatabilityAllometric
Agate.Configuration.AssimilationBinary
```

### Normalization and process compilation

```@docs
Agate.Processes.NamedProcess
Agate.Processes.NormalizedModelDefinition
Agate.Processes.ParameterRequirementIdentity
Agate.Processes.ParameterRequirement
Agate.Processes.ParameterBinding
Agate.Processes.ParameterApplicability
Agate.Processes.driver_identities
Agate.Processes.parameter_requirements
Agate.Processes.parameter_bindings
Agate.Processes.parameter_name
Agate.Processes.resolve_parameter_applicability
Agate.Processes.normalize_model
Agate.Processes.process_rate
Agate.Compilation.model_contributions
Agate.Compilation.group_contributions
Agate.Compilation.compile_tendency
Agate.Compilation.compile_tendencies
Agate.Compilation.compile_model_tendencies
```

## Named families, recipes, and replay

Named model families add stable code identity and durable recipe replay around the same definition-driven process compiler. `ProcessModelRecipe` is the `agate.model_recipe.v3` scientific representation for component/process families. `ModelManifest` records the resolved execution state.

```@docs
Agate.Construction.ProcessModelRecipe
Agate.Construction.ModelManifest
Agate.Construction.construct_factory
Agate.Construction.construct_factory_plus_manifest
Agate.Construction.capture_process_model_recipe
Agate.Construction.replay_factory
Agate.Construction.resolve_construction_scalar_type
Agate.Construction.recipe_family
Agate.Construction.recipe_factory
Agate.Construction.encode_recipe
Agate.Construction.decode_recipe
Agate.Construction.export_recipe
Agate.Construction.import_recipe
Agate.Factories.AbstractBGCFactory
Agate.Factories.parameter_definitions
Agate.Factories.parameter_directory
Agate.Factories.parameter_spec
Agate.Factories.default_components
Agate.Factories.default_processes
Agate.Factories.default_community
Agate.Configuration.PFTSpecification
Agate.Configuration.build_plankton_community
```

`ModelRecipe` represents the v2 recipe format; component/process models use `ProcessModelRecipe` v3.

```@docs
Agate.Construction.ModelRecipe
```

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

pal.rows
pal.columns
pal.matrix
```
