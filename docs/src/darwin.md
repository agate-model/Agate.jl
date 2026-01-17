### [Agate.jl DARWIN model](@id DARWIN)

Agate's simplified DARWIN-like model is constructed via `DARWIN.construct`.

The constructor follows the same high-level pattern as `NiPiZD.construct`, but
exposes an additional hook for overriding selected biogeochemical tracer dynamics
by key.

```@docs
Agate.DARWIN.construct
```
