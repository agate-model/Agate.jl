# Equation-based dynamics API

Agate builds biogeochemical kernels from *equations*: construction-time symbolic objects that

1. assemble a plain Julia `Expr` for each tracer tendency, and
2. record which **parameters** and **interaction matrices** are required.

At runtime (CPU or GPU), the model executes a normal `Expr` containing only arithmetic on numeric arrays and scalars. The symbolic objects never enter kernels.

## Parameter references: `PV.<name>`

Dynamics authors write equations using **parameter placeholders** passed in as the first argument `PV`.

Dynamics builders accept a first argument `PV` (provided by `construct`). Reference parameters through `PV`:

- Scalars: `PV.detritus_remineralization`
- Per-size-class vectors: `PV.maximum_predation_rate[i]`
- Interaction matrices: `PV.palatability_matrix[j, i]`, `PV.assimilation_matrix[j, i]`

These `PV.<name>` bindings are small placeholders (`ParamVar`) that record requirements when indexed.

During `construct`, Agate builds `PV` from the active parameter registry. So when you add new parameters, you typically **extend the registry** (with `ParamSpec`s) and keep writing equations using `PV.<name>` — no extra declaration step.

## Missing / `nothing` policy

Agate no longer encodes missing behaviour in the *equation syntax*.
Instead, each `ParamSpec` in the parameter registry declares how missing values are handled during CPU resolution:

- `missing_policy = :fail` — missing/`nothing` is an error
- `missing_policy = :zero_warn` — replace with `0`/`false` and warn
- `missing_policy = :zero_silent` — replace with `0`/`false` silently

This keeps the equation authoring surface clean and GPU-safe.

## `sum_over`: symbolic sum builder

`sum_over` expands a list of terms into a plain `Expr` sum at construction time.

```julia
function mortality_sum(PV, plankton_syms)
    return sum_over(plankton_syms) do sym, i
        PV.linear_mortality[i] * sym
    end
end
```

This is a construction-time helper only: loops are unrolled into the final expression so kernels remain GPU-friendly.

## `Equation`

All plankton and biogeochemical dynamics builders must return an `Equation`.

- `Equation` cannot be constructed from a raw `Expr`.
- Use library building blocks (or parameter identifiers and `sum_over`) to assemble expressions.

Example: a detritivorous heterotroph growth tendency

```julia
using Agate.Equations: Equation
using Agate.Library.Mortality: linear_loss
using Agate.Library.Predation: grazing_loss

@inline monod(D, k) = D / (D + k)

function heterotroph_growth(PV, plankton_syms::AbstractVector{Symbol}, plankton_sym::Symbol, plankton_idx::Int)
    uptake = PV.maximum_detritus_uptake_rate[plankton_idx] *
             monod(:D, PV.detritus_half_saturation[plankton_idx]) *
             plankton_sym

    grazing = grazing_loss(PV, plankton_sym, plankton_idx, plankton_syms)
    mort = linear_loss(PV, plankton_sym, plankton_idx)

    return Equation(uptake - grazing - mort)
end
```
