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

bgc_type = construct(NiPiZDFactory())
bgc = bgc_type()
```

If you need a different float type (for example, GPU execution), specify `FT` explicitly:

```julia
bgc_type_f32 = construct(NiPiZDFactory(); FT=Float32)
```

## Structural defaults vs parameter defaults

Factories provide *structural* defaults (community size structure) and default dynamics functions.
All *parameter defaults* are sourced from the model's parameter registry (single source of truth).

- Structural defaults: `Agate.Models.default_parameter_args(factory)`
- Dynamics defaults: `Agate.Models.default_plankton_dynamics(factory)`, `Agate.Models.default_biogeochem_dynamics(factory)`
- Parameter defaults: `Agate.Parameters.parameter_registry(factory)`

## Overriding community and parameters

Use `default_parameter_args` to apply readable keyword overrides.

```julia
using Agate
using Agate.Constructor: construct, default_parameter_args, update_plankton_args, update_dynamics
using Agate.Models: NiPiZDFactory
using Agate.Library.Allometry: AllometricParam, PowerLaw
using Oceananigans.Units: day

factory = NiPiZDFactory()

# 1) Customize community structure (sizes/diameters).
community = Agate.Models.default_parameter_args(factory)
community = update_plankton_args(community, :Z; n=1, diameters=[60.0])
community = update_plankton_args(community, :P; n=3, diameters=(1.5, 20.0, :log_splitting))

# 2) Optionally swap in alternative dynamics builders.
plankton_dynamics = Agate.Models.default_plankton_dynamics(factory)
biogeochem_dynamics = Agate.Models.default_biogeochem_dynamics(factory)

# 3) Parameter overrides (registry keys as keywords).
parameter_args = default_parameter_args(factory;
    community,
    params=(
        detritus_remineralization = 0.18 / day,
        maximum_growth_rate = (P = AllometricParam(PowerLaw(); prefactor=3.0 / day, exponent=-0.15),),
    ),
)

bgc_type = construct(factory;
    plankton_dynamics,
    biogeochem_dynamics,
    parameter_args=parameter_args,
    sinking_tracers=(D = 2.0 / day,),
)

bgc = bgc_type()
```

## API reference

```@docs
Agate.Constructor.construct
Agate.Constructor.default_parameter_args
Agate.Constructor.patch
Agate.Constructor.update_plankton_args
Agate.Constructor.update_dynamics
```
