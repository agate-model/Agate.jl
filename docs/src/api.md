# API reference

## Direct model construction

`ModelDefinition` is the scientific model-definition container. Direct construction realizes intrinsic component size structure, resolves process-owned parameter requirements and required drivers, and compiles runtime tracer equations during setup.

```@docs
Agate.Construction.construct
Agate.Processes.ModelDefinition
```

### Components

```@docs
Agate.Configuration.Population
Agate.Configuration.Pool
Agate.Configuration.currency
Agate.Configuration.states
Agate.Configuration.state_currency
Agate.Configuration.size_structure
Agate.Configuration.ComponentLayout
Agate.Configuration.realize_components
Agate.Configuration.realize_component_groups
Agate.Configuration.component_classes
Agate.Configuration.component_state_tracers
Agate.Configuration.state_tracers
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
Agate.Processes.Consumption
Agate.Processes.Grazing
Agate.Processes.Mortality
Agate.Processes.ProductRouting
Agate.Processes.FixedStoichiometry
Agate.Processes.Remineralization
```

### Parameter definitions

```@docs
Agate.Parameters.ParameterSpec
Agate.Parameters.ParameterProvision
Agate.Parameters.ParameterDefinition
Agate.Parameters.DefaultProvider
Agate.Parameters.ConstDefault
Agate.Parameters.DerivedDefault
Agate.Parameters.derive_default
Agate.Parameters.NoDefault
Agate.Parameters.FillDefault
Agate.Parameters.DiameterIndexedVectorDefault
Agate.Parameters.DiameterIndexedMaterialization
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
Agate.Compilation.compile_model_tendencies
```

## Named families, recipes, and replay

Named model families add stable code identity and durable recipe replay around the same definition-driven process compiler. `ProcessModelRecipe` is the `agate.model_recipe.v0.4` scientific representation for component/process families. `ModelManifest` records the resolved execution state.

```@docs
Agate.Construction.ProcessModelRecipe
Agate.Construction.ModelManifest
Agate.Construction.construct_plus_manifest
Agate.Construction.capture_process_model_recipe
Agate.Construction.replay_family
Agate.Construction.resolve_construction_scalar_type
Agate.Construction.family_id
Agate.Construction.registered_family
Agate.Construction.encode_recipe
Agate.Construction.decode_recipe
Agate.Construction.export_recipe
Agate.Construction.import_recipe
Agate.ModelFamilies.AbstractModelFamily
Agate.Parameters.parameter_definitions
Agate.Parameters.parameter_directory
Agate.Parameters.parameter_spec
Agate.ModelFamilies.default_components
Agate.ModelFamilies.default_processes
Agate.Configuration.PFTSpecification
Agate.Configuration.build_plankton_community
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
