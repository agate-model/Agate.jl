### [Agate.jl-NiPiZD model](@id NiPiZD)

NiPiZD uses two phytoplankton SizeClasses (`2` and `10` μm) and two zooplankton SizeClasses
(`20` and `100` μm) by default. The default predator:prey diameter ratio is `10`, and living-prey
grazing uses `PreferentialGrazing(switching_exponent=1)`, so each consumer shares one maximum
ingestion capacity across its prey community.

Selected default biological parameters are:

| Parameter | Default |
|---|---:|
| maximum growth allometry | `1 / day`, exponent `-0.15` |
| maximum grazing allometry | `21.9 / day`, exponent `-0.16` |
| linear mortality | `0.02 / day` |
| quadratic zooplankton mortality | `0` |
| assimilation efficiency | `0.7` |
| grazing half-saturation | `16 / 106 mmol N m^-3` |

```@docs
Agate.Models.NiPiZD.construct
Agate.Models.NiPiZD.construct_plus_recipe
Agate.Models.NiPiZD.construct_from_recipe
```
