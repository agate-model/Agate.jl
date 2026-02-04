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

##### Example: Z cannibalism toggle

To *allow* Z–Z interactions (cannibalism), include `:Z` on the prey axis. Then you can enable/disable
cannibalism by writing the `:Z` prey columns to zero or to a nonzero value.

```julia
using Agate

# Build a community where Z is both a consumer and a potential prey.
roles = (consumers=(:Z,), prey=(:P, :Z))

# Disable cannibalism by zeroing the prey columns whose group is :Z.
pal_no_cannibalism = (ctx) -> begin
    nC = length(ctx.consumer_indices)
    nP = length(ctx.prey_indices)
    P = ones(ctx.FT, nC, nP)

    prey_groups = ctx.group_symbols[ctx.prey_indices]
    z_cols = findall(==(:Z), prey_groups)

    P[:, z_cols] .= zero(ctx.FT)  # Z does not eat Z
    return P
end

bgc = NiPiZD.construct(; roles, palatability_matrix=pal_no_cannibalism)

# Alternatively, if cannibalism should be *impossible* (and you want smaller matrices),
# omit :Z from the prey axis entirely:
bgc_small = NiPiZD.construct(; roles=(consumers=(:Z,), prey=(:P,)))
```

Roles define which plankton classes participate as **consumers** and **prey**. By default,
NiPiZD uses consumers `(:Z,)` and prey `(:P,)`, but you can override this via the
`roles=(consumers=..., prey=...)` keyword (group symbols, indices, or boolean masks).

Interaction matrices are stored and used as **rectangular consumer-by-prey matrices**.
If you do not provide overrides, Agate derives default matrices from the model parameters
for the chosen roles.

```@docs
Agate.NiPiZD.construct
```
