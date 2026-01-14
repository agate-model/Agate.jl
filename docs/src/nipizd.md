### [Agate.jl NiPiZD model](@id NiPiZD)

Agate's NiPiZD model is configured through a **factory** (`NiPiZDFactory`) and constructed via
`Agate.Constructor.construct`.

At a high level:

1. Choose a factory (here, `NiPiZDFactory()`).
2. Optionally override components (dynamics functions), community structure (sizes/diameters),
   or parameter values (via `default_parameter_args(...; params=(...))`).
3. Call `construct(factory; ...)` to obtain a concrete biogeochemistry type.
4. Instantiate the type and pass it to Oceananigans/OceanBioME.

See **Examples → Box model factories** for a worked, end-to-end box model demonstration of
component swapping and parameter overrides.

```@docs
Agate.Models.NiPiZDFactory
Agate.Constructor.construct
```
