### [Agate.jl NiPiZD model](@id NiPiZD)

Agate's NiPiZD model is constructed via `NiPiZD.construct`.

The constructor is intentionally small and explicit:

- **Community structure**: choose `n_phyto`, `n_zoo`, and diameter specifications
- **Dynamics**: optionally swap any of the four default dynamics builders
- **Parameters**: override named parameters via `parameters=(; ...)`
- **Interactions**: optionally provide `palatability_matrix` and `assimilation_matrix` overrides (matrices or provider functions)

The returned value is an Oceananigans/OceanBioME-compatible biogeochemistry instance.

```@docs
Agate.NiPiZD.construct
```
