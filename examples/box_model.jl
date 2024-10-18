# This example shows how to create a Biogeochemistry (BGC) model from parameter and tracer
# definitions and run a box model (0D) simulation using OceanBioME and Oceananigans.

using Agate
using Plots

# ==================================================
# Define BGC model (NPZD)
# ==================================================

model_path = joinpath("NPZD", "model.jl")
include(model_path)
model = NPZD()

# ==================================================
# Run box model
# ==================================================

init_conditions = (N=7.0, P=0.01, Z=0.05, D=0.0)
timeseries = run_box_model(model, init_conditions)

# ==================================================
# Plotting
# ==================================================

p = plot(timeseries.P; label="P")
plot!(p, timeseries.Z; label="Z")
plot!(p, timeseries.D; label="D")
savefig(p, "NPZD_box_oceananigans.png")
