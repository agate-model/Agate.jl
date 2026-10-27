### [Agate.jl-Onyx model](@id NiPiZD)

Onyx uses two phytoplankton SizeClasses (`2` and `10` μm) and two zooplankton SizeClasses
(`20` and `100` μm) by default. The default predator:prey diameter ratio is `10`, and living-prey
grazing uses `PreferentialGrazing(switching_exponent=1)`, so each consumer shares one maximum
ingestion capacity across its prey community.

```@docs
Onyx.construct
Onyx.construct_plus_recipe
```
