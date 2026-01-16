# Constructor API

This page summarizes Agate's **factory-based constructor surface**.

Agate models are configured via **factories** (for example `NiPiZDFactory()` and `DarwinFactory()`) and then
compiled into a concrete biogeochemistry **instance** using `construct`.

## Core workflow

`construct` returns an *instance* (already `Adapt.jl` compatible):

```julia
using Agate

bgc = construct(NiPiZDFactory())  # CPU instance by default
```

If you need a different float type (for example, for GPU execution), specify `FT` explicitly:

```julia
bgc_f32 = construct(NiPiZDFactory(); FT=Float32)
```

To adapt the constructed instance to a GPU backend, either:

- adapt yourself (explicit and transparent), or
- construct directly on a GPU architecture via `arch=GPU()`.

```julia
using Adapt
using CUDA
using Oceananigans.Architectures: GPU

bgc = construct(NiPiZDFactory())
bgc_gpu = Adapt.adapt(CuArray, bgc)

# or (recommended when using Oceananigans/OceanBioME)
bgc_gpu2 = construct(NiPiZDFactory(); FT=Float32, arch=GPU())
```

## Introspection helpers

For interactive work (REPL/notebooks), Agate provides small helpers to discover
the model surface:

```julia
names = tracer_names(bgc)                 # Vector{Symbol}
pars  = parameter_names(bgc)              # Vector{Symbol}
```

## Structural defaults vs parameter defaults

Factories provide *structural* defaults (community size structure) and default dynamics functions.
All *parameter defaults* are sourced from the model's parameter registry (single source of truth).

- Structural defaults: `default_community(factory)`
- Dynamics defaults: `Agate.Models.default_plankton_dynamics(factory)`, `Agate.Models.default_biogeochem_dynamics(factory)`
- Parameter defaults: `parameter_registry(factory)`

## Overriding community and parameters

Community structure is overridden via `update_community`, while parameter values are overridden by updating the **registry**.

```julia
using Agate
using Agate.Library.Allometry: AllometricParam, PowerLaw
using Oceananigans.Units: day

factory = NiPiZDFactory()

# 1) Customize community structure (sizes/diameters).
community = default_community(factory)
community = update_community(community, :Z; n=1, diameters=[60.0])
community = update_community(community, :P; n=3, diameters=(1.5, 20.0, :log_splitting))

# 2) Override parameters by updating the registry.
registry = parameter_registry(factory)
registry = update_registry(registry;
    detritus_remineralization = 0.18 / day,
    maximum_growth_rate = (P = AllometricParam(PowerLaw(); prefactor=3.0 / day, exponent=-0.15),),
)

# 3) Compile the model instance.
bgc = construct(factory;
    community,
    registry,
    sinking_tracers=(D = 2.0 / day,),
)
```

## Interactions

Interaction matrices (and other interaction-related overrides) flow through the same registry mechanism.
Provide them via the `interactions` keyword as either:

- a `NamedTuple` of overrides, or
- a function `(ctx) -> NamedTuple` (useful when matrices depend on community structure).

Unknown keys should error (typo protection).

## API reference

```@docs
Agate.construct
Agate.default_community
Agate.tracer_names
Agate.parameter_names
Agate.required_parameters
Agate.update_community
Agate.extend_community
Agate.parameter_registry
Agate.update_registry
Agate.extend_registry
Agate.update_dynamics
Agate.extend_dynamics
Agate.define_tracer_functions
Agate.Constructor.patch
```
