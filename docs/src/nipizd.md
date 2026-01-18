### [Agate.jl NiPiZD model](@id NiPiZD)

Agate's NiPiZD model is constructed via `NiPiZD.construct`.

The constructor is intentionally small and explicit:

- **Community structure**: choose `n_phyto`, `n_zoo`, and diameter specifications
- **Dynamics**: optionally swap any of the four default dynamics builders
- **Parameters**: override named parameters via `parameters=(; ...)`
- **Interactions**: optionally provide palatability and assimilation matrices

  - simplest: pass concrete matrices via `palatability_matrix=` and/or `assimilation_matrix=`
  - beginner-friendly provider (no `ctx ->`): pass a function
    `f(diameters, group_symbols) -> matrix`
  - advanced: pass a context provider `f(ctx) -> matrix`

The returned value is an Oceananigans/OceanBioME-compatible biogeochemistry instance.

```@docs
Agate.NiPiZD.construct
```
