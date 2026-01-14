### [Agate.jl DARWIN model](@id DARWIN)

Agate's DARWIN-like model is configured through a **factory** (`DarwinFactory`) and constructed via
`Agate.Constructor.construct`.

The factory-based API is consistent across models:

1. Choose a factory (here, `DarwinFactory()`).
2. Optionally override community structure (sizes/diameters), tracer dynamics (`plankton_dynamics` and
   `biogeochem_dynamics`), or parameter values (via `default_parameter_args(...; params=(...))`).
3. Call `construct(factory; ...)` to obtain a concrete biogeochemistry type.
4. Instantiate the type and use it with Oceananigans/OceanBioME.

The **paper/GPU** scripts show an end-to-end example of running a DARWIN configuration on the GPU.

```@docs
Agate.Models.DarwinFactory
Agate.Constructor.construct
```
