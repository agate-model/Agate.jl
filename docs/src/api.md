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

`ProcessModelRecipe` uses the `agate.model_recipe.v3` schema for component/process model
families. It records canonical components, named processes, semantic parameter bindings,
subgroup realization, and authored overrides; runtime precision, host fields, parameter
materialization, topology maps, and compiled equations are derived during realization.
`ModelRecipe` remains the `agate.model_recipe.v2` representation used by model families
that still realize through the existing factory dynamics path. `ModelManifest` records
the resolved execution state for either construction path.

```@docs
Agate.Construction.ModelRecipe
Agate.Construction.ProcessModelRecipe
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
Agate.Construction.capture_process_model_recipe
Agate.Construction.replay_factory
Agate.Construction.resolve_construction_scalar_type
Agate.Construction.recipe_family
Agate.Construction.recipe_factory
Agate.Factories.AbstractBGCFactory
Agate.Factories.parameter_definitions
Agate.Factories.parameter_directory
Agate.Factories.parameter_spec
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
Agate.Factories.default_components
Agate.Factories.default_processes
Agate.Factories.default_community
Agate.Factories.default_plankton_dynamics
Agate.Factories.default_biogeochem_dynamics
Agate.Configuration.PalatabilityAllometric
Agate.Configuration.AssimilationBinary
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
Agate.Configuration.realize_component_groups
```

### Processes and model definitions

```@docs
Agate.Processes.AbstractProcess
Agate.Processes.AbstractFormulation
Agate.Processes.AbstractFactor
Agate.Processes.Growth
Agate.Processes.Light
Agate.Processes.NutrientResponse
Agate.Processes.Temperature
Agate.Processes.factors
Agate.Processes.Grazing
Agate.Processes.Consumption
Agate.Processes.Mortality
Agate.Processes.ProductRouting
Agate.Processes.Remineralization
Agate.Processes.ModelDefinition
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
```

### Process contribution compilation

```@docs
Agate.Compilation.model_contributions
Agate.Compilation.group_contributions
Agate.Compilation.compile_tendency
Agate.Compilation.compile_tendencies
Agate.Compilation.compile_model_tendencies
```

### Plankton communities and interactions

```@docs
Agate.Configuration.PFTSpecification
Agate.Configuration.build_plankton_community
Agate.Configuration.PalatabilityAllometric
Agate.Configuration.AssimilationBinary
```
