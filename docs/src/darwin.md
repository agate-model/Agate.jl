### [Agate.jl DARWIN model](@id DARWIN)

Agate's simplified DARWIN-like model is configured through a **factory** (`DarwinFactory`) and constructed via
`construct`.

At a high level:

1. Choose a factory (here, `DarwinFactory()`).
2. Optionally override:
   - community structure via `default_community` + `update_community`
   - parameter values via `parameter_registry` + `update_registry`
   - dynamics builders via `plankton_dynamics` / `biogeochem_dynamics` keywords to `construct`
3. Call `construct(factory; ...)` to obtain a concrete biogeochemistry **instance**.
4. Pass the instance to Oceananigans/OceanBioME (and adapt to GPU if desired).

```@docs
Agate.DarwinFactory
Agate.construct
Agate.default_community
Agate.parameter_registry
Agate.update_community
Agate.update_registry
```
