# Editing interaction blocks by group

When you want to quickly change **who eats whom** at the group level (e.g. toggle zooplankton cannibalism), it is convenient to work with a small group-by-group matrix.

Agate provides two small helpers:

  - `roles_from_groups(...)` to define consumer/prey membership by group `Symbol`s.
  - `interaction_blocks(...)` to build an editable group-block matrix (with explicit group ordering).

## Example: Z eats P, optional Z cannibalism

```julia
using Agate
using Agate.Library.Light
using Agate.Configuration: roles_from_groups, interaction_blocks, set_block!, forbid_link!

# Allow Z to be both a consumer and (optionally) prey.
roles = roles_from_groups(; consumers=:Z, prey=(:P, :Z))

# Create an editable group-by-group block matrix in the same group order as `roles`.
pal = interaction_blocks(roles; init=0)

# Z eats P strongly.
set_block!(pal; consumer_group=:Z, prey_group=:P, value=1.0)

# Start with cannibalism enabled (then toggle it off).
set_block!(pal; consumer_group=:Z, prey_group=:Z, value=0.25)
forbid_link!(pal; consumer_group=:Z, prey_group=:Z)

# Next, build a model instance using these role axes and overrides.
```

`NiPiZD.construct` keeps its surface small and does not expose role customization. To use custom role axes, call the model-agnostic constructor:

```julia
using Agate
using OceanBioME: BoxModelGrid
using Agate.Configuration: roles_from_groups, interaction_blocks, set_block!, forbid_link!

factory = Agate.Models.NiPiZD.NiPiZDFactory()
base = Agate.Factories.default_community(factory)
community = Agate.Construction.build_plankton_community(
    base;
    n=(Z=2, P=2),
    diameters=(Z=(20, 100, :linear_splitting), P=(2, 10, :log_splitting)),
)

bgc = Agate.Construction.construct_factory(
    factory;
    grid=BoxModelGrid(),
    community=community,
    interaction_roles=roles,
    interaction_overrides=(palatability_matrix=pal,),
    auxiliary_fields=(:PAR,),
)
```

Notes:

  - `pal` is normalized during construction into the canonical rectangular consumer-by-prey matrix.
  - The block helpers are host-side utilities; they do not run inside kernels.

## Scaling links

```julia
using Agate.Configuration: scale_block!

scale_block!(pal; consumer_group=:Z, prey_group=:P, factor=0.5)
```

## Passing a full group-block matrix

If you have a full group-by-group matrix over *all* groups (in the explicit order of the `community` `NamedTuple` used for construction), wrap it explicitly:

```julia
using Agate.Configuration: GroupBlockMatrix

B = [0.0 1.0; 0.0 0.0]
bgc = NiPiZD.construct(; palatability_matrix=GroupBlockMatrix(B))
```
