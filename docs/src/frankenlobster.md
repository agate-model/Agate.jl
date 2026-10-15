### [Agate.jl-FrankenLOBSTER model](@id FrankenLOBSTER)

FrankenLOBSTER generates an size structured plankton community that can then be used inside OceanBioME's `LOBSTER` nutrient and detritus model. By default, it has two phytoplankton (`P_1`, `P_2`), two zooplankton (`Z_1`, `Z_2`), and one heterotrophic
bacterioplankton (`H_1`) type. Phytoplankton use nitrate and ammonium,
bacterioplankton consume dissolved organic matter, and zooplankton graze on
living phytoplankton and bacterioplankton.

The [FrankenLOBSTER with OceanBioME example](generated/frankenlobster.md)
shows how to compose the plankton component with `LOBSTER`, run a short
water-column simulation, and capture and replay a model recipe.

```@docs
Agate.Models.FrankenLOBSTER.construct
Agate.Models.FrankenLOBSTER.construct_plus_recipe
```
