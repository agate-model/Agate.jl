### [Agate.jl NiPiZD model](@id NiPiZD)

Agate's NiPiZD model is configured through a **factory** (`NiPiZDFactory`) and constructed via
`Agate.Models.construct`.

At a high level:

1. Choose a factory (here, `NiPiZDFactory()`).
2. Optionally override components (dynamics functions), community structure (`plankton_args`),
   or parameters (PFT and biogeochemistry specifications).
3. Call `construct(factory; ...)` to obtain a concrete biogeochemistry type.
4. Instantiate the type and pass it to Oceananigans/OceanBioME.

See **Examples → Box model factories** for a worked, end-to-end box model demonstration of
component swapping and parameter overrides.

```@docs
Agate.Models.NiPiZDFactory
Agate.Models.construct
```
