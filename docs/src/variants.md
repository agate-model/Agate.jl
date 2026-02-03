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

## Adding a new variant (developer)

The recommended pattern is:

 1. Create a file under `src/Models/<FAMILY>/Variants/`.
 2. Implement a builder `(; kwargs...) -> VariantSpec`.
 3. Call `register_variant(id, builder)` during module initialization (usually in the variant file).

A minimal sketch:

```julia
using Agate.Models: ModelId, VariantSpec, register_variant
using Agate.Models.DARWIN: DarwinFactory
using Agate.Interface:
    default_plankton_dynamics, default_biogeochem_dynamics, default_community

function citation2026_B_spec(; n_phyto=2, n_zoo=2)
    id = ModelId(:DARWIN, :citation2026, :B)
    factory = DarwinFactory()

    return VariantSpec(
        id,
        factory,
        default_plankton_dynamics(factory),
        default_biogeochem_dynamics(factory),
        default_community(factory),
        (;),
        nothing,
    )
end

register_variant(ModelId(:DARWIN, :citation2026, :B), citation2026_B_spec)
```

Keep variant defaults minimal: store only the parameters/interaction_overrides that differ from the family defaults.
