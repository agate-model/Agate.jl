# Constructor API

This page summarises Agate's **factory-based constructor surface**.

Agate models are configured via **factories** (for example `NiPiZDFactory()` and `DarwinFactory()`) and then
compiled into a concrete biogeochemistry type using `Agate.Constructor.construct`.

## Core workflow

`construct` returns a **type**, which you then instantiate:

```julia
using Agate
using Agate.Constructor: construct
using Agate.Models: NiPiZDFactory

bgc_type = construct(NiPiZDFactory()) # defaults to FT = Float64
bgc = bgc_type()
```

If you need a different float type (for example, GPU execution), specify `FT` explicitly:

```julia
bgc_type_f32 = construct(NiPiZDFactory(); FT=Float32)
```

## Overriding configuration

Factories expose defaults for plankton community structure, dynamics, and non-plankton parameters via
`Agate.Models.default_*` functions. You can pull those defaults, then override them and pass the
modified objects back into `construct`.

Agate provides two small helpers for ergonomic, non-mutating overrides:

- `patch`: override fields on immutable configuration objects (`NamedTuple`s and parameter containers).
- `update_group`: convenience helper for overriding a single plankton group entry inside `plankton_args`.

Example:

```julia
using Agate
using Agate.Models: NiPiZDFactory
using Agate.Constructor: construct, patch, update_group

using Oceananigans.Units: day

factory = NiPiZDFactory()

plankton_args = Agate.Models.default_plankton_args(factory, Float64)
biogeochem_args = Agate.Models.default_biogeochem_args(factory, Float64)

#show defaults _args

# 1) Override community structure.
plankton_args = update_group(plankton_args, :Z; n=1, diameters=[60.0])
plankton_args = update_group(plankton_args, :P; n=3, diameters=(1.5, 20.0, :log_splitting))

# 2) Override a PFT parameter (creates a new PFTSpecification).
phyto_pft_fast = patch(plankton_args.P.pft; maximum_growth_rate_a=3.0/day)
plankton_args = update_group(plankton_args, :P; pft=phyto_pft_fast)

# 3) Override a biogeochemistry specification value.
biogeochem_args = patch(biogeochem_args; detritus_remineralization=0.18/day)

# 4) Construct a concrete biogeochemistry type.
bgc_type = construct(factory; plankton_args=plankton_args, biogeochem_args=biogeochem_args)
bgc = bgc_type()
```

## API reference

```@docs
Agate.Constructor.construct
Agate.Constructor.patch
Agate.Constructor.update_group
Agate.Utils.Specifications.PFTSpecification
Agate.Utils.Specifications.BiogeochemistrySpecification
Agate.Utils.Specifications.ModelSpecification
```
