# Equation-based dynamics API

Agate builds biogeochemical kernels from *equations*: construction-time symbolic objects that

1. assemble a plain Julia `Expr` for each tracer tendency, and
2. record which **parameters**, **interaction matrices**, and **biogeochemical scalars** are required.

At runtime (CPU or GPU), the model still executes a normal `Expr` containing only arithmetic on numeric arrays and scalars. The symbolic objects never enter the kernels.

## Public names (set C)

These names are the user-facing, construction-time API for dynamics authors:

- `Equation(...)` — wrapper returned by all dynamics builders.
- `group_param(:key)[i]` — a parameter used by a group **in its own dynamics**.
- `community_param(:key)[i]` — a parameter that may be inactive for some groups (missing or `nothing` is allowed).
- `interaction_matrix(:key)[j, i]` — an interaction matrix entry.
- `bgc_param(:key)` — a biogeochemical *scalar* coming from the `BiogeochemistrySpecification`.
- `Σ(items) do sym, idx ... end` — construction-time symbolic sum builder (unrolled into an `Expr`).

## Requirements and missing vs `nothing`

Agate enforces a fixed missing/`nothing` policy during construction:

- **Explicit `nothing`** means *inactive*: the constructor fills **zeros** for the affected indices.
- **Missing key** (typo protection):
  - if a key is required by a group’s own dynamics (referenced via `group_param`) and is missing from that group’s PFT spec → **error**.
  - if a key is only referenced as a community-optional parameter (via `community_param`) → no error; the constructor fills **zeros**.

This policy is not configurable.

## What each reference means

### `group_param(:key)`

Use when the parameter must be defined (or explicitly set to `nothing`) for the group whose equation you are building.

Example:

```julia
mort = linear_loss(P, i)  # internally uses group_param(:linear_mortality)[i]
```

### `community_param(:key)`

Use when the parameter can be undefined/inactive for some groups, and missing should not be an error.

Typical use: predator parameters referenced in prey loss sums.

```julia
loss = grazing_loss(prey, prey_i, plankton_syms)
# internally uses community_param(:maximum_predation_rate)[pred_j], etc.
```

### `bgc_param(:key)`

Use for scalar fields (not per-group vectors) supplied in `biogeochem_args::BiogeochemistrySpecification`.

```julia
export_frac = bgc_param(:mortality_export_fraction)
remin       = bgc_param(:detritus_remineralization) * :D
```

`bgc_param` is distinct from `community_param` because it is a **scalar** (shared across the whole model),
while `community_param` is a **vector parameter container** indexed by plankton size-class.

### `interaction_matrix(:key)`

Use for interaction matrices (e.g. palatability and assimilation efficiency).

```julia
p = interaction_matrix(:palatability_matrix)[pred_j, prey_i]
```

If a required matrix is not provided explicitly, Agate builds default matrices using library builders (based on PFT traits), with optional user overrides.

## `Σ`: symbolic sum builder

`Σ` expands a list of terms into a plain `Expr` sum at construction time.

```julia
loss_sum = Σ(plankton_syms) do sym, i
    community_param(:linear_mortality)[i] * sym
end
```

This is a construction-time helper only: loops are unrolled into the final expression so kernels remain GPU-friendly.

## `Equation`

All plankton and biogeochemical dynamics builders must return an `Equation`.

- `Equation` cannot be constructed from a raw `Expr`.
- Use library building blocks (or `group_param` / `community_param` / `interaction_matrix` / `bgc_param` / `Σ`) to assemble expressions.

Example: phytoplankton default

```julia
using Agate.Library.Equations: Equation
using Agate.Library.Mortality: linear_loss
using Agate.Library.Predation: grazing_loss
using Agate.Library.Photosynthesis: growth_single_nutrient

function phytoplankton_default(plankton_syms, P::Symbol, i::Int)
    growth  = growth_single_nutrient(:N, P, :PAR, i)
    grazing = grazing_loss(P, i, plankton_syms)
    mort    = linear_loss(P, i)
    return Equation(growth - grazing - mort)
end
```

