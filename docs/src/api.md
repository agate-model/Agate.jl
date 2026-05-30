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

```@docs
Agate.Models.construct_with_manifest
Agate.Models.export_manifest
Agate.Construction.PFTSpecification
Agate.Factories.AbstractBGCFactory
Agate.Equations.CompiledEquation
Agate.Construction.define_tracer_functions
```

## Plankton communities and size structure

```@docs
Agate.Configuration.build_plankton_community
```
