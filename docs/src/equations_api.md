# Parameters and interaction matrices

Agate model constructors are explicit: you choose community structure and dynamics, then optionally override
**parameters** and **interactions**.

This page documents the parameter and interaction surface that is shared by the built-in models.

## Parameter metadata: `parameter_directory`

Each model factory defines a **parameter directory** via `Agate.Utils.parameter_directory(factory)`.
The directory is a collection of `ParameterSpec` entries describing:

  - the parameter key (a `Symbol`)
  - the expected **shape** (`:scalar`, `:vector`, `:matrix`)
  - the value kind (`:real` or `:bool`)
  - for matrices, optional role-aware `axes` (e.g. `(:consumer, :prey)`)

The constructor uses this metadata to validate overrides early (typo protection + shape checks).

## Overriding parameters

All model constructors accept a `parameters=(; ...)` NamedTuple.

```julia
bgc = NiPiZD.construct(; parameters=(; detritus_remineralization=0.18 / day,))
```

Validation behaviour:

  - Unknown keys throw immediately.
  - Vectors must have length `n_total` (one value per plankton class).
  - Matrices must match the declared `ParameterSpec` shape (details below).

### Interaction traits

NiPiZD and DARWIN also expose a small set of *interaction traits* (vectors) that are used to derive
default interaction matrices:

  - `can_eat`, `can_be_eaten` (Bool)
  - `optimum_predator_prey_ratio`, `specificity`, `protection` (Real)
  - `assimilation_efficiency` (Real)

These are validated like any other vector parameter.

## Interaction matrices

Models expose two interaction matrices:

  - `palatability_matrix`
  - `assimilation_matrix`

These are **role-aware** predator-by-prey matrices with axes `(:consumer, :prey)`.

### Accepted override forms

For either matrix key, you may pass:

  - a full square matrix of size `(n_total, n_total)`
  - a rectangular axis-sized matrix of size `(n_consumer, n_prey)`
  - an axis-local group-block matrix of size `(n_consumer_groups, n_prey_groups)`
  - a group-block matrix over *all* groups, wrapped as `Agate.Utils.GroupBlockMatrix(B)`
  - a provider function `(ctx) -> matrix` returning any of the above

Role-aware rectangular matrices are the preferred override form because they are explicit and small.

### The construction context: `InteractionContext`

Provider functions receive an `Agate.Utils.InteractionContext` with precomputed axes:

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

bgc = NiPiZD.construct(; palatability_matrix=pal)
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

bgc = NiPiZD.construct(; n_phyto, n_zoo, parameters=(; specificity=fill(0.15f0, n_total),))
# palatability_matrix is regenerated from the updated specificity
```

### Access at runtime

To keep kernels GPU-friendly, Agate stores role-aware matrices canonically as rectangular
`(n_consumer, n_prey)` arrays under:

```julia
bgc.parameters.interactions.palatability_matrix
bgc.parameters.interactions.assimilation_matrix
```

For convenience, Agate also exposes **square views** on the top-level parameter NamedTuple:

```julia
bgc.parameters.palatability_matrix
bgc.parameters.assimilation_matrix
```

These square views are indexed by *global* plankton indices and return zero outside the consumer/prey axes.
