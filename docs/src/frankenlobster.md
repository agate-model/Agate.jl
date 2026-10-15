### [Agate.jl-FrankenLOBSTER model](@id FrankenLOBSTER)

FrankenLOBSTER constructs a size-structured plankton community for composition with
OceanBioME's `LOBSTER` nutrient and detritus model. Agate represents the living community,
while `LOBSTER` supplies nitrate (`NO₃`), ammonium (`NH₄`), dissolved organic matter (`DOM`), 
and small and large particulate organic matter (`sPOM`, `bPOM`).
Photosynthetically active radiation (`PAR`) is supplied by the light-attenuation model.

The default community has two phytoplankton SizeClasses (`P_1`, `P_2`, diameters 0.6 and
1.2 μm), two zooplankton (`Z_1`, `Z_2`, 6 and 12 μm), and one heterotrophic
bacterioplankton SizeClass (`H_1`, 0.6 μm). Phytoplankton take up nitrate and ammonium,
with ammonium inhibiting nitrate uptake; bacteria consume `DOM`; and zooplankton graze
on living phytoplankton and bacteria. The coupling returns nutrient uptake and waste
fluxes to OceanBioME's native pools.

## Construction and settings

`construct` returns an `Agate.Integrations.NPDPlankton` component, ready to pass to
`LOBSTER`. The `size_structure` keyword controls the phytoplankton, zooplankton and
bacterioplankton PFTs. Use `parameters` for biological rate and interaction overrides;
use `settings` for model-level conversion factors such as `chlorophyll_ratio` (default
`1.31`) and `carbon_ratio` (default `6.56`). The resolved settings are available
as `plankton.settings`. The package requires OceanBioME `0.19`.

```julia
using Agate
using OceanBioME.Models.NutrientsPlanktonDetritusModels: LOBSTER
using Oceananigans

grid = RectilinearGrid(CPU(); size=(1, 1, 4), extent=(1, 1, 40))
plankton = Agate.Models.FrankenLOBSTER.construct(;
    grid,
    parameters=(maximum_growth_rate=(P_1=1e-5,),),
    settings=(chlorophyll_ratio=1.5,),
)
bgc = LOBSTER(grid; plankton)
```

The [FrankenLOBSTER with OceanBioME example](generated/frankenlobster.md) runs a short
water-column simulation and demonstrates recipe capture and replay.

## Recipes

`construct_plus_recipe` returns the plankton component and a versioned `ModelRecipe`.
The `agate.model_recipe.v0.2` document captures the family definition version, size
structure, parameters, settings and sinking configuration. Export a recipe as JSON and
reconstruct it with the execution grid of your choice:

```julia
plankton, recipe = Agate.Models.FrankenLOBSTER.construct_plus_recipe(; grid)
Agate.Construction.export_recipe("frankenlobster.json", recipe)
replayed = Agate.Models.FrankenLOBSTER.construct(
    Agate.Construction.import_recipe("frankenlobster.json"); grid
)
```

## API

```@docs
Agate.Models.FrankenLOBSTER.construct
Agate.Models.FrankenLOBSTER.construct_plus_recipe
```
