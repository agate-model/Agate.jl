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

Agate provides a few small helpers for ergonomic, non-mutating overrides:

- `update_plankton_args`: update a single plankton group entry inside `plankton_args` (and optionally patch its `pft`).
- `update_biogeochem_args`: update biogeochemistry specification values.
- `update_dynamics`: update a dynamics `NamedTuple` with key validation.

Example:

```julia
using Agate
using Agate.Models: NiPiZDFactory
using Agate.Constructor: construct, update_plankton_args, update_biogeochem_args
using Agate.Library.Allometry: AllometricParam, PowerLaw

using Oceananigans.Units: day

factory = NiPiZDFactory()

plankton_args  = Agate.Models.default_plankton_args(factory)
biogeochem_args = Agate.Models.default_biogeochem_args(factory)

# 1) Override community structure.
plankton_args = update_plankton_args(plankton_args, :Z; n=1, diameters=[60.0])
plankton_args = update_plankton_args(plankton_args, :P; n=3, diameters=(1.5, 20.0, :log_splitting))

# 2) Override a PFT parameter (one step).
plankton_args = update_plankton_args(plankton_args, :P;
    maximum_growth_rate = AllometricParam(PowerLaw(); prefactor=3.0/day, exponent=-0.15),
)

# 3) Override a biogeochemistry specification value.
biogeochem_args = update_biogeochem_args(biogeochem_args; detritus_remineralization=0.18/day)

# 4) Construct a concrete biogeochemistry type.
bgc_type = construct(factory; plankton_args=plankton_args, biogeochem_args=biogeochem_args)
bgc = bgc_type()
```

## API reference

```@docs
Agate.Constructor.construct
Agate.Constructor.update_plankton_args
Agate.Constructor.update_biogeochem_args
Agate.Constructor.update_dynamics
Agate.Constructor.patch
Agate.Utils.Specifications.PFTSpecification
Agate.Utils.Specifications.BiogeochemistrySpecification
Agate.Utils.Specifications.ModelSpecification
```
