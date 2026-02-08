# Parameters and interaction matrices

Agate model constructors are explicit: you choose community structure and dynamics, then optionally override
**parameters** and **interactions**.

This page documents the parameter and interaction surface that is shared by the built-in models.

## Parameter metadata: `parameter_directory`

Each model factory defines a **parameter directory** via `Agate.Factories.parameter_directory(factory)`.
The directory is a collection of `ParameterSpec` entries describing:

  - the parameter key (a `Symbol`)
  - the expected **shape** (`:scalar`, `:vector`, `:matrix`)
  - for matrices, optional role-aware `axes` (e.g. `(:consumer, :prey)`)

The constructor uses this metadata to validate overrides early (typo protection + shape checks).

For developers: `parameter_directory(factory)` is derived from `parameter_definitions(factory)` by default, so specs and defaults can be declared in one place.

## Overriding parameters

All model constructors accept a `parameters=(; ...)` NamedTuple.

```julia
bgc = Agate.Models.NiPiZD.construct(; parameters=(; detritus_remineralization=0.18 / day,))
```

Validation behaviour:

  - Unknown keys throw immediately.
  - Vectors must have length `n_total` (one value per plankton class).
  - Matrices must match the declared `ParameterSpec` shape (details below).

### Interaction traits

NiPiZD and DARWIN expose a small set of *interaction traits* (vectors) that are used to derive
default interaction matrices over the active consumer/prey role axes:

  - `optimum_predator_prey_ratio`, `specificity`, `protection` (Real)
  - `assimilation_efficiency` (Real)

Role axes are defined via the `interaction_roles` argument to `Agate.Construction.construct_factory` (typically using group symbols).

These traits are validated like any other vector parameter.

When you override any of these traits *without* explicitly overriding the corresponding interaction matrix, Agate recomputes the derived matrix during construction so the resolved parameter set remains consistent.

## Interaction matrices

Models expose two interaction matrices:

  - `palatability_matrix`
  - `assimilation_matrix`

These are **role-aware** predator-by-prey matrices with axes `(:consumer, :prey)`.

### Accepted override forms

For either matrix key, you may pass:

  - a rectangular axis-sized matrix of size `(n_consumer, n_prey)`
  - a provider function `(ctx) -> matrix` returning one of the above

Full-square `(n_total, n_total)` matrices are not accepted as overrides for role-aware parameters.
If you want to override based on a full matrix, slice it to the active axes yourself.

Role-aware rectangular matrices are the preferred override form because they are explicit and small.

### The construction context: `CommunityContext`

Provider functions receive an `Agate.Configuration.CommunityContext` with precomputed axes:

  - `ctx.consumer_indices` / `ctx.prey_indices` (global indices)
  - `ctx.group_symbols`, `ctx.diameters`, `ctx.n_total`
  - `ctx.FT` (the floating-point type inferred from the grid)

Example: override palatability as a rectangular matrix:

```julia
pal = (ctx) -> begin
    nc = length(ctx.consumer_indices)
    np = length(ctx.prey_indices)

    M = zeros(ctx.FT, nc, np)
    # fill M in consumer-by-prey order
    return M
end

bgc = Agate.Models.NiPiZD.construct(; palatability_matrix=pal)
```

### Derived matrices (trait-driven)

During `construct`, Agate resolves defaults, then applies overrides.
If a matrix is **not** explicitly overridden, but at least one of its declared *trait dependencies*
*is* explicitly overridden, Agate recomputes the matrix.

This makes it practical to tune interaction behaviour via a small, readable set of vectors.

```julia
n_phyto = 4
n_zoo = 2
n_total = n_phyto + n_zoo

bgc = Agate.Models.NiPiZD.construct(; n_phyto, n_zoo, parameters=(; specificity=fill(0.15f0, n_total),))
# palatability_matrix is regenerated from the updated specificity
```

### Access at runtime

Agate stores role-aware matrices canonically as rectangular `(n_consumer, n_prey)` arrays:

```julia
bgc.parameters.palatability_matrix
bgc.parameters.assimilation_matrix
```

The same matrices (plus axis mappings) are available in the internal helper container:

```julia
ints = bgc.parameters.interactions
ints.palatability
ints.assimilation
ints.consumer_global      # axis → global
ints.prey_global          # axis → global
ints.global_to_consumer   # global → axis (0 means "not on axis")
ints.global_to_prey       # global → axis (0 means "not on axis")
```
