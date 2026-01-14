# Equation-based dynamics API

Agate builds biogeochemical kernels from *equations*: construction-time symbolic objects that

1. assemble a plain Julia `Expr` for each tracer tendency, and
2. record which **parameters** and **interaction matrices** are required.

At runtime (CPU or GPU), the model executes a normal `Expr` containing only arithmetic on numeric arrays and scalars. The symbolic objects never enter kernels.

## Parameter references: bare identifiers

Dynamics authors write equations using *bare parameter identifiers*:

- Scalars: `detritus_remineralization`
- Per-group vectors: `maximum_predation_rate[i]`
- Interaction matrices: `palatability_matrix[j, i]`

These identifiers are small placeholders (`ParamVar`) that record requirements when indexed. Each model module automatically declares these placeholders from its parameter registry.

If you are authoring dynamics outside a model module and you need **new** parameter names, you can declare placeholders locally:

```julia
using Agate.Library.Equations: @paramvars
@paramvars maximum_detritus_uptake_rate detritus_half_saturation
```

## Missing / `nothing` policy

Agate no longer encodes missing behaviour in the *equation syntax*.
Instead, each `ParamSpec` in the parameter registry declares how missing values are handled during CPU resolution:

- `scope = :fail` — missing/`nothing` is an error
- `scope = :zero_warn` — replace with `0`/`false` and warn
- `scope = :zero_silent` — replace with `0`/`false` silently

This keeps the equation authoring surface clean and GPU-safe.

## `Σ`: symbolic sum builder

`Σ` expands a list of terms into a plain `Expr` sum at construction time.

```julia
loss_sum = Σ(plankton_syms) do sym, i
    linear_mortality[i] * sym
end
```

This is a construction-time helper only: loops are unrolled into the final expression so kernels remain GPU-friendly.

## `Equation`

All plankton and biogeochemical dynamics builders must return an `Equation`.

- `Equation` cannot be constructed from a raw `Expr`.
- Use library building blocks (or parameter identifiers and `Σ`) to assemble expressions.

Example: a detritivorous heterotroph growth tendency

```julia
using Agate.Library.Equations: Equation
using Agate.Library.Mortality: linear_loss
using Agate.Library.Predation: grazing_loss

@inline monod(D, k) = D / (D + k)

function heterotroph_growth(plankton_syms::AbstractVector{Symbol}, plankton_sym::Symbol, plankton_idx::Int)
    uptake = maximum_detritus_uptake_rate[plankton_idx] *
             monod(:D, detritus_half_saturation[plankton_idx]) *
             plankton_sym

    grazing = grazing_loss(plankton_sym, plankton_idx, plankton_syms)
    mort = linear_loss(plankton_sym, plankton_idx)

    return Equation(uptake - grazing - mort)
end
```
