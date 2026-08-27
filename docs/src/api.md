# API reference

## Direct model construction

`ModelDefinition` is the scientific model-definition container. Direct construction realizes intrinsic component size structure, resolves process-owned parameter slots and required drivers, and compiles runtime tracer equations during setup.

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
Agate.Configuration.population_state
Agate.Configuration.size_structure
Agate.Configuration.component_classes
Agate.Configuration.component_state_tracers
Agate.Configuration.state_tracers
Agate.Configuration.state_tracer
```

### Processes and factors

Formulations are authored as concrete scientific objects, for example `Light(Smith(); driver=:PAR)`
and `Consumption(HeterotrophicConsumption(); ...)`. Numerical scientific parameters belong to
the model parameter system rather than the formulation object; for example
`Nutrients(FrankTNorm(); ...)` uses the Frank t-norm's declared `sharpness` parameter slot.
`FrankTNorm()` names the formulation; `Agate.Library.Nutrients.frank_tnorm` is the numerical
kernel. Parameterized nodes bind their formulation-local slots directly to model-level parameter
names with `bindings=(...)`. Omitted slots bind by the same name; a `Symbol` explicitly renames
or shares one parameter, while a one-level qualifier map handles repeated slots such as
source-specific remineralization. Formulation and factor authoring is method-based rather than
registry-based. Custom process topologies can extend Agate through the narrow
`Processes.normalize_process_facts` and `Compilation.process_fluxes` hooks, keeping custom
topology in the same normalization and construction pipeline as built-in processes. Durable recipes
identify a registered family and its version rather than serializing process or formulation objects.

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
Agate.Processes.Mortality
Agate.Processes.Products
Agate.Processes.FixedStoichiometry
Agate.Processes.Remineralization
Agate.Processes.authored_parameter_bindings
```

### Parameter definitions

The keyed parameter block separates runtime process parameters from construction-only inputs.
Scientific slots and realized process applicability determine `Parameter` vector or matrix storage
automatically, so runtime parameters never restate axes. `MetaParameter` values exist only during
construction to feed `DerivedDefault` calculations; shaped meta-parameters declare their setup
domain with `axes=`. Scientific slot-to-parameter relationships are authored beside the process
or factor through `bindings=`.

```@docs
Agate.Parameters.Parameter
Agate.Parameters.MetaParameter
Agate.Parameters.DefaultProvider
Agate.Parameters.DerivedDefault
Agate.Parameters.derive_default
Agate.Parameters.NoDefault
Agate.Parameters.DiameterIndexedVectorDefault
Agate.Configuration.PalatabilityAllometric
Agate.Configuration.AssimilationBinary
```

### Normalization and process compilation

```@docs
Agate.Processes.NamedProcess
Agate.Processes.NormalizedModelDefinition
Agate.Processes.ParameterBinding
Agate.Processes.driver_identities
Agate.Processes.parameter_bindings
Agate.Processes.normalize_model
Agate.Compilation.compile_model_tendencies
```

## Named families, recipes, and replay

Named model families add stable code identity and durable recipe replay around the same definition-driven process compiler. `ProcessModelRecipe` is the `agate.model_recipe.v1` family/version/realization document: it records the registered family, exact `definition_version`, canonical population/size realization, parameter overrides, sinking choices, and bottom state. The loaded family supplies the canonical component/process definition on replay. `ModelManifest` records the resolved execution state.

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
Agate.ModelFamilies.definition_version
Agate.Parameters.parameter_definitions
Agate.ModelFamilies.default_components
Agate.ModelFamilies.default_processes
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
pal = interaction_matrix(bgc, :palatability_matrix)

pal.rows
pal.columns
pal.matrix
```
