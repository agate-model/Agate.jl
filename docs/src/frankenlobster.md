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

`construct` returns an `Agate.Integrations.NPDPlankton` component ready for
`LOBSTER`. Use `size_structure` to choose the phytoplankton, zooplankton and
bacterioplankton PFTs; `parameters` overrides biological rates and interactions,
while `settings` overrides model-level ratios. The resolved defaults are
`chlorophyll_ratio=1.31` (nonnegative) and `carbon_ratio=6.56` (positive).
They can be inspected with `plankton.settings` or
`Agate.Introspection.model_settings(plankton)`.

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

PFT names are canonicalized in name order and each explicit diameter list in
ascending order. Realized SizeClasses are named `<PFT>_1`, `<PFT>_2`, etc.,
so positional parameter overrides follow that realized order. **Named** overrides
such as `(maximum_growth_rate=(P_1=1e-5,),)` identify their target explicitly.

## Coupling and tracers

Agate owns the living plankton (`P_1`, `P_2`, `Z_1`, `Z_2`, `H_1` by default).
`LOBSTER` supplies native nutrient (`NO₃`, `NH₄`) and detritus (`DOM`, `sPOM`,
`bPOM`) components. Agate's configured nutrient tracers must be a subset of
native nutrient tracers; unconfigured native nutrients have zero Agate plankton
uptake. Consumed detritus must be present in the native detritus component.
The adapter validates these at composition. By default, the exchange channels map
`:solid` to `:solid_waste`, `:dissolved` to `:dissolved_waste`, and `:inorganic`
to `:inorganic_waste` when those runtime pools exist, or to `nothing` otherwise.
A custom low-level adapter specifies all three channel keys, in any order, with
`nothing` for an unused channel. NPD model families supply both conversion ratios
in their resolved model settings; standalone runtimes can supply them through the
`settings` keyword. Both are visible through `Agate.Introspection.model_settings`.

`Agate.Introspection.tracer_names(plankton)` lists the compiled Agate runtime's
tracers, including external nutrient and exchange-pool names. In contrast,
`required_biogeochemical_tracers(plankton)` identifies only Agate-owned living
tracers required from the plankton component. The fully coupled tracer set is
available from `required_biogeochemical_tracers(bgc.underlying_biogeochemistry)`.

Temperature (`T`) is supplied to the coupled model as a **prognostic tracer**
and read by Agate as a physical dependency. Initialize it in the Oceananigans
model, for example with `set!(model; T=20.0, ...)` as in the
[FrankenLOBSTER with OceanBioME example](generated/frankenlobster.md).
Photosynthetically active radiation (`PAR`) is instead an auxiliary field from
the chosen light-attenuation model.

Waste/exchange pools are flux destinations, evaluated through OceanBioME's signed
coupling hooks; their concentrations are read as zero when assembling the Agate
runtime state. The current nutrient, detritus and waste exchange hooks evaluate
the biological tendencies with `t=0`, using the instantaneous tracer and
auxiliary fields. The compiled runtime remains responsible for the biological
equations.

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
