### [Agate.jl DARWIN model](@id DARWIN)

Agate's simplified DARWIN-like model is constructed via `DARWIN.construct`.

The high-level pattern matches `NiPiZD.construct`, with additional elemental cycling tracers.

- **Community structure**: choose `n_phyto`, `n_zoo`, and diameter specifications
- **Dynamics**: optionally swap plankton dynamics and override selected biogeochemical tracer dynamics by key
- **Parameters**: override named parameters via `parameters=(; ...)`
- **Interactions**: optionally override palatability and assimilation matrices

DARWIN shares the same interaction-matrix surface as NiPiZD: rectangular `(n_consumer, n_prey)` overrides,
provider functions `(ctx) -> matrix`, and trait-driven derived matrices.

```@docs
Agate.DARWIN.construct
```
