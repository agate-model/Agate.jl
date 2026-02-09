### [Agate.jl NiPiZD model](@id NiPiZD)

Agate's NiPiZD model is constructed via `Agate.Models.NiPiZD.construct`.

The constructor is intentionally small and explicit:

  - **Community structure**: choose `n_phyto`, `n_zoo`, and diameter specifications
  - **Dynamics**: optionally swap any of the four default dynamics builders
  - **Parameters**: override named parameters via `parameters=(; ...)`
  - **Interactions**: optionally override palatability and assimilation matrices

#### Interaction overrides

Interaction matrices are role-aware predator-by-prey matrices.
The preferred override form is a rectangular `(n_consumer, n_prey)` matrix.

You can provide interaction overrides as concrete matrices (`palatability_matrix=...`, `assimilation_matrix=...`).

Provider functions / callables are not supported in user overrides. If you need matrices derived from traits or other parameters, define a `Variant` / `Factory` default that produces concrete rectangular matrices during construction.


```@docs
Agate.Models.NiPiZD.construct
```
