### [Agate.jl NiPiZD model](@id NiPiZD)

Agate's NiPiZD model is constructed via `NiPiZD.construct`.

The constructor is intentionally small and explicit:

- **Community structure**: choose `n_phyto`, `n_zoo`, and diameter specifications
- **Dynamics**: optionally swap any of the four default dynamics builders
- **Parameters**: override named parameters via `parameters=(; ...)`
- **Interactions**: optionally override palatability and assimilation matrices

#### Interaction overrides

Interaction matrices are role-aware predator-by-prey matrices.
The preferred override form is a rectangular `(n_consumer, n_prey)` matrix.

You can provide interaction overrides as:

- concrete matrices (`palatability_matrix=...`, `assimilation_matrix=...`)
- provider functions `(ctx) -> matrix`, evaluated once during construction

NiPiZD also exposes a small set of interaction *traits* (vectors) that are used to derive default
interaction matrices. If you override one of these traits (and do not explicitly override the
corresponding matrix), Agate recomputes the matrix during construction.

```@docs
Agate.NiPiZD.construct
```
