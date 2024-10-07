# This example shows how to create a Biogeochemistry (BGC) model from parameter and tracer
# definitions and run a box model (0D) simulation using OceanBioME and Oceananigans.

using Agate
using Plots

# ==================================================
# Define BGC model (NPZD)
# ==================================================

# NPZD `parameters` and `tracers` are defined here
include("NPZD/model_definition.jl")

BGC = create_bgc_struct(:BGC, parameters)
add_bgc_methods(
    BGC,
    tracers,
    auxiliary_fields=aux_field_vars,
    helper_functions="NPZD/functions.jl"
)
bgc_model = BGC()

# ==================================================
# Run box model
# ==================================================

init_conditions = (N = 7.0, P = 0.01, Z = 0.05, D=0.0)
timeseries = run_box_model(bgc_model, init_conditions)

# ==================================================
# Plotting
# ==================================================

p = plot(timeseries.P, label="P")
plot!(p, timeseries.Z, label="Z")
plot!(p, timeseries.D, label="D")
savefig(p, "NPZD_box_oceananigans.png")
