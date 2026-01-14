"""Agate.ParamVars

Canonical namespace for construction-time parameter placeholders.

Equation authoring should reference parameters through this module, e.g.

```julia
const PV = Agate.ParamVars

uptake = PV.maximum_growth_rate[i] * P
```

The placeholders themselves are declared programmatically from the active
parameter registry during `construct`.

Notes
-----
- This module is intentionally *empty* by default.
- Placeholders are `ParamVar{:name}()` values from `Agate.Library.Equations`.
- The declaration step should happen once per construction run, based on the
  merged registry.
"""
module ParamVars

end # module
