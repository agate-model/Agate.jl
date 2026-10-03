# # [FrankenLOBSTER with OceanBioME] (@id frankenlobster_example)

# FrankenLOBSTER is Agate's richer P/Z/H model-family example. Agate owns the living
# phytoplankton, zooplankton, and heterotrophic-bacterioplankton ecology; OceanBioME's
# `LOBSTER` supplies the surrounding nutrient and detritus model.

using Agate
using OceanBioME
using OceanBioME.Models.NutrientsPlanktonDetritusModels: LOBSTER
using Oceananigans
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Fields: ConstantField
using Oceananigans.Units: minutes, hours

# ## Construct and compose the ecosystem

# The public constructor returns an OceanBioME-compatible plankton component backed by a
# compiled Agate runtime. Model-family parameters are overridden through `parameters`, while
# model-level conversion/diagnostic properties use `settings`.

grid = RectilinearGrid(CPU(); size=(1, 1, 4), extent=(1, 1, 40))
plankton = Agate.Models.FrankenLOBSTER.construct(; grid)

# Constant PAR keeps this short example focused on biological composition rather than physical
# forcing.
light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100.0))
bgc = LOBSTER(grid; plankton, light_attenuation)

# The coupled tracer set combines Agate's living community with OceanBioME-owned nutrient,
# detritus, and temperature state.
required_biogeochemical_tracers(bgc.underlying_biogeochemistry)

# ## Run a short water-column experiment

model = NonhydrostaticModel(grid; biogeochemistry=bgc)
set!(
    model;
    NO₃=7.0,
    NH₄=0.1,
    DOM=0.1,
    sPOM=0.01,
    bPOM=0.01,
    T=20.0,
    P_1=0.01,
    P_2=0.01,
    Z_1=0.02,
    Z_2=0.02,
    H_1=0.01,
)

simulation = Simulation(model; Δt=5minutes, stop_time=12hours)
run!(simulation)

# The example is intentionally short: it demonstrates end-to-end construction and execution,
# rather than serving as a calibrated ecological experiment.
living_tracers = (:P_1, :P_2, :Z_1, :Z_2, :H_1)
final_living = NamedTuple{living_tracers}(
    Tuple(getproperty(model.tracers, tracer)[1, 1, 1] for tracer in living_tracers)
)
final_living

# ## Capture and replay the scientific realization

# Recipe replay uses the same `construct` function. Julia dispatch distinguishes a new authored
# construction from a replay based on the `ModelRecipe` argument.
_, recipe = Agate.Models.FrankenLOBSTER.construct_plus_recipe(; grid)
replayed = Agate.Models.FrankenLOBSTER.construct(recipe; grid)

required_biogeochemical_tracers(replayed)
