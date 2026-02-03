# Editing interaction blocks by group

When you want to quickly change **who eats whom** at the group level (e.g. toggle zooplankton cannibalism), it is convenient to work with a small group-by-group matrix.

Agate provides two small helpers:

- `roles_from_groups(...)` to define consumer/prey membership by group `Symbol`s.
- `interaction_blocks(...)` to build an editable group-block matrix (with explicit group ordering).

## Example: Z eats P, optional Z cannibalism

```julia
using Agate
using Agate.Library.Light
using Agate.Models.NiPiZD: NiPiZDFactory
using Agate.Utils: roles_from_groups, interaction_blocks, set_block!, forbid_link!

# Allow Z to be both a consumer and (optionally) prey.
roles = roles_from_groups(consumers = :Z, prey = (:P, :Z))

# Create an editable group-by-group block matrix in the same group order as `roles`.
pal = interaction_blocks(NiPiZDFactory(); roles)

# Z eats P strongly.
set_block!(pal; consumer_group = :Z, prey_group = :P, value = 1.0)

# Start with cannibalism enabled (then toggle it off).
set_block!(pal; consumer_group = :Z, prey_group = :Z, value = 0.25)
forbid_link!(pal; consumer_group = :Z, prey_group = :Z)

bgc = NiPiZD.construct(; roles, palatability_matrix = pal)
```

Notes:

- `pal` is normalized during construction into the canonical rectangular consumer-by-prey matrix.
- The block helpers are host-side utilities; they do not run inside kernels.

## Scaling links

```julia
using Agate.Utils: scale_block!

scale_block!(pal; consumer_group = :Z, prey_group = :P, factor = 0.5)
```

## Passing a full group-block matrix

If you have a full group-by-group matrix over *all* groups (in the order returned by `required_groups(factory)`), wrap it explicitly:

```julia
using Agate.Utils: GroupBlockMatrix

B = [0.0 1.0; 0.0 0.0]
bgc = NiPiZD.construct(; palatability_matrix = GroupBlockMatrix(B))
```