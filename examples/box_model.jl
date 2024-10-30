# This example shows how to create a Biogeochemistry (BGC) model from parameter and tracer
# definitions and run a box model (0D) simulation using OceanBioME and Oceananigans.

using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Plots

# ==================================================
# Define BGC model (NPZD with cyclical PAR)
# ==================================================

include(joinpath("NPZD", "tracers.jl"))
bgc_model = Biogeochemistry(NPZD(); light_attenuation=FunctionPAR(; grid=BoxModelGrid()))

# ==================================================
# Run box model
# NOTE: this could be an Oceananigans model instead
# ==================================================

full_model = BoxModel(; biogeochemistry=bgc_model)
init_conditions = (N=7.0, P=0.01, Z=0.05, D=0.0)

timeseries = run_simulation(full_model, init_conditions)

# ==================================================
# Plotting
# ==================================================

p = plot(timeseries.P; label="P")
plot!(p, timeseries.Z; label="Z")
plot!(p, timeseries.D; label="D")
savefig(p, "NPZD_box_oceananigans.png")
