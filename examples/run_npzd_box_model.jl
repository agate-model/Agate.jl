# This example shows how to create a model from parameter and tracer definitions and run
# a box model (0D) simulation.

using Agate
using Plots

# NPZD `parameters` and `tracers` are defined here
include("NPZD/model_definition.jl")

# create NPZD model
NPZD = create_bgc_struct(:NPZD, parameters)
add_bgc_methods(
    NPZD,
    tracers,
    auxiliary_fields=aux_field_vars,
    helper_functions="NPZD/functions.jl"
)
npzd_model = NPZD()

# set initial conditions and run box model simulation
init_conditions = (N = 7.0, P = 0.01, Z = 0.05, D=0.0)
timeseries = run_box_model(npzd_model, init_conditions)

# plot results
p = plot(timeseries.P, label="P")
plot!(p, timeseries.Z, label="Z")
plot!(p, timeseries.D, label="D")
savefig(p, "NPZD.png")
