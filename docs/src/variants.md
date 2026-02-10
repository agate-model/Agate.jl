# Variants

Agate supports **model variants**: named, reproducible recipes for constructing a model.

A variant is meant to capture choices that you typically want to keep fixed for a manuscript or a long-running experiment, such as:

  - which plankton / biogeochemical **dynamics** functions are used,
  - default **community structure** (groups, size classes, diameters, PFT traits),
  - default **parameter overrides** and **interaction overrides**.

## Identifying a variant

Variants are identified by a `ModelId(family, citation, tag)`:

  - `family`: the model family, e.g. `:NiPiZD` or `:DARWIN`.
  - `citation`: a stable key for the manuscript or configuration lineage, e.g. `:citation2026`.
  - `tag`: a stable label for a concrete recipe within that citation, e.g. `:A`, `:B`, `:submission`, `:accepted`.

## Listing variants

```julia
using Agate
using Agate.Models: list_variants

list_variants()
list_variants(; family=:DARWIN)
```

## Constructing from a variant

Use `variant(...)` to obtain a `VariantSpec`, then construct it.

```julia
using Agate
using Agate.Models: ModelId, variant, construct

id = ModelId(:DARWIN, :citation2026, :A)

# The builder may accept kwargs such as n_phyto/n_zoo to set the community.
spec = variant(id; n_phyto=2, n_zoo=2)

bgc = construct(spec; parameters=(;), interaction_overrides=nothing, grid=BoxModelGrid())
```

### Runtime overrides

  - `parameters=...` is merged on top of the spec's default `spec.parameters`.
  - `interaction_overrides=...` is merged on top of the spec's default `spec.interaction_overrides`.

This makes variants a good place to define *defaults*, while keeping per-run changes explicit.

## Deriving interaction matrices

User-facing interaction overrides are **data-only**: you must pass explicit rectangular matrices.
If you want interaction matrices to be *derived* from traits or other parameters, do that work in a
Variant/Factory layer and ensure the factory receives concrete rectangular matrices during construction.

The built-in DARWIN factory derives `palatability_matrix` and `assimilation_matrix` from trait vectors.
If you want to change the derivation algorithm, do it in a custom Variant/Factory by overriding
`matrix_definitions(::MyFactory)` to swap the *strategy object* used for a given matrix.

```julia
using Agate
using Agate.Models: ModelId, variant, construct

spec = variant(ModelId(:DARWIN, :citation2026, :A); n_phyto=8, n_zoo=4)
bgc = construct(spec; grid=BoxModelGrid())

size(bgc.parameters.palatability_matrix)   # (n_consumer, n_prey)
size(bgc.parameters.assimilation_matrix)  # (n_consumer, n_prey)
```

### Python/MATLAB-style workflow

1. Choose a built-in variant or create your own.
2. Construct the model.
3. If you need to change interactions at runtime, pass **explicit rectangular matrices**.

```julia
spec = variant(:DARWIN, :citation2026, :A)
bgc0 = construct(spec; grid=BoxModelGrid())

# Override interactions with explicit matrices (no function-valued overrides).
n_cons, n_prey = size(bgc0.parameters.palatability_matrix)
M = fill(0.1, n_cons, n_prey)
bgc1 = construct(spec; grid=BoxModelGrid(), interaction_overrides=(; palatability_matrix=M,))
```

## Adding a new variant (developer)

The recommended pattern is:

 1. Create a file under `src/Models/<FAMILY>/Variants/`.
 2. Implement a builder `(; kwargs...) -> VariantSpec`.
 3. Call `register_variant(id, builder)` during module initialization (usually in the variant file).

A minimal sketch:

```julia
using Agate.Models: ModelId, VariantSpec, register_variant
using Agate.Models.DARWIN: DarwinFactory
using Agate.Factories:
    default_plankton_dynamics, default_biogeochem_dynamics, default_community

function citation2026_B_spec(; n_phyto=2, n_zoo=2)
    id = ModelId(:DARWIN, :citation2026, :B)
    factory = DarwinFactory()

    interaction_roles = (consumers=(:Z,), prey=(:P,))
    auxiliary_fields = (:PAR,)

    return VariantSpec(
        id,
        factory,
        default_plankton_dynamics(factory),
        default_biogeochem_dynamics(factory),
        default_community(factory),
        interaction_roles,
        auxiliary_fields,
        (;),
        nothing,
    )
end

register_variant(ModelId(:DARWIN, :citation2026, :B), citation2026_B_spec)
```

Keep variant defaults minimal: store only the parameters/interaction_overrides that differ from the family defaults.
