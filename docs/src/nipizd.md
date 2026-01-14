### [Agate.jl NiPiZD model](@id NiPiZD)

Agate's NiPiZD model is configured through a **factory** (`NiPiZDFactory`) and constructed via
`construct`.

At a high level:

1. Choose a factory (here, `NiPiZDFactory()`).
2. Optionally override:
   - community structure (sizes/diameters) via `default_community` + `update_community`
   - parameter values via `parameter_registry` + `update_registry`
   - dynamics builders via `plankton_dynamics` / `biogeochem_dynamics` keywords to `construct`
3. Call `construct(factory; ...)` to obtain a concrete biogeochemistry **instance**.
4. Pass the instance to Oceananigans/OceanBioME (and adapt to GPU if desired).

See **Examples → Box model factories** for a worked end-to-end box model demonstration of
component swapping and parameter overrides.

```@docs
Agate.NiPiZDFactory
Agate.construct
Agate.default_community
Agate.parameter_registry
Agate.update_community
Agate.update_registry
```
