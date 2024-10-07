using Agate
using Plots

include("NPZD_definition.jl")

NPZD = create_bgc_struct(:NPZD, parameters)
add_bgc_methods(
    NPZD,
    tracers,
    auxiliary_fields=aux_field_vars,
    helper_functions="NPZD_functions.jl"
)
npzd_model = NPZD()

init_conditions = (N = 7.0, P = 0.01, Z = 0.05, D=0.0)
timeseries = run_boxmodel(npzd_model, init_conditions)

p = plot(timeseries.P, label="P")
plot!(p, timeseries.Z, label="Z")
plot!(p, timeseries.D, label="D")
savefig(p, "NPZD.png")
