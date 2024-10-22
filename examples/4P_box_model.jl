# This example shows how to create a Biogeochemistry (BGC) model from parameter and tracer
# definitions and run a box model (0D) simulation using OceanBioME and Oceananigans.

using Agate
using Plots
using Oceananigans.Units

const year = years = 365day

using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields

# ==================================================
# Define BGC model (NPZD)
# ==================================================

include("N2P2ZD/model.jl")
model = N2P2ZD()

# ==================================================
# Run box model
# ==================================================

print(required_biogeochemical_tracers(model))

#init_conditions = (P1=0.1, P2=0.1, Z1=0.05, Z2=0.05, N=7.0, D=1)
init_conditions = (N=7.0, Z2=0.05, D=0.0, P1=0.1, P2=0.1, Z1=0.05)

timeseries = run_box_model(
    model, 
    Δt=10minutes,
    stop_time=10years,   
    init_conditions)

# ==================================================
# Plotting
# ==================================================

p1 = plot(timeseries.P1; label="P1", title="Phytoplankton 1")
p2 = plot(timeseries.P2; label="P2", title="Phytoplankton 2")
p3 = plot(timeseries.Z1; label="Z1", title="Zooplankton 1")
p4 = plot(timeseries.Z2; label="Z2", title="Zooplankton 2")
p5 = plot(timeseries.D; label="D", title="Detritus")
p6 = plot(timeseries.N; label="N", title="Nitrogen")

p = plot(p1, p2, p3, p4, p5, p6; layout=(2, 3), size=(900, 600))

savefig(p, "N2P2ZD_box_oceananigans.png")