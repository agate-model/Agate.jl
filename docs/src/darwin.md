### [Agate.jl DARWIN model](@id DARWIN)

Agate's simplified DARWIN-like model is constructed via `Agate.Models.DARWIN.construct`.

The high-level pattern matches `Agate.Models.NiPiZD.construct`, with additional elemental cycling tracers.

  - **Community structure**: choose `n_phyto`, `n_zoo`, and diameter specifications
  - **Dynamics**: optionally swap plankton dynamics and override selected biogeochemical tracer dynamics by key
  - **Parameters**: override named parameters via `parameters=(; ...)`
  - **Interactions**: optionally override palatability and assimilation matrices

DARWIN shares the same interaction-matrix surface as NiPiZD: rectangular `(n_consumer, n_prey)` overrides and trait-driven derived matrices.

Provider functions / callables are not supported in user overrides. If you need matrices derived from traits or other parameters, define a `Variant` / `Factory` default that produces concrete rectangular matrices during construction.

See also: the Variants section for how to define a Variant/Factory that changes interaction-matrix defaults or derivation strategies via `matrix_definitions`.

```@docs
Agate.Models.DARWIN.construct
```
