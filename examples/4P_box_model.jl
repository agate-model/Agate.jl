# This example shows how to create a Biogeochemistry (BGC) model from parameter and tracer
# definitions and run a box model (0D) simulation using OceanBioME and Oceananigans.

using Agate
using Plots

# ==================================================
# Define BGC model (NPZD)
# ==================================================

include("N2P2ZD/model.jl")
model = N2P2ZD()

# ==================================================
# Run box model
# ==================================================
init_conditions = (P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)

timeseries = run_box_model(model, init_conditions)

# ==================================================
# Plotting
# ==================================================

p = plot(timeseries.P1; label="P1")
plot!(p, timeseries.Z1; label="Z1")
plot!(p, timeseries.D; label="D")
savefig(p, "N2P2ZD_box_oceananigans.png")
