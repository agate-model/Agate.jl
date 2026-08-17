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

### Recipes and realizations

```@docs
Agate.Construction.ModelRecipe
Agate.Construction.ModelRealization
Agate.Construction.encode_recipe
Agate.Construction.decode_recipe
Agate.Construction.export_recipe
Agate.Construction.import_recipe
```

## Plankton communities and size structure

```@docs
Agate.Configuration.build_plankton_community
```
