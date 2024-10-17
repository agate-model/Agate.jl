# This example shows how to integrate a Biogeochemistry (BGC) box model with DifferentialEquations.
using Agate.Light

using DifferentialEquations
using Plots

using Oceananigans.Units
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields

const year = years = 365day

# ==================================================
# Define BGC model (NPZD)
# ==================================================

include("N2P2ZD/model.jl")
model = N2P2ZD()

# ==================================================
# DifferentialEquations
# ==================================================

function model_ODEs(du, u, p, t)
    model = N2P2ZD(p...)

    # NOTE: in more complex examples there might be other auxiliary fields that should be
    # calculated here and passed to the function below
    PAR = cyclical_PAR(t)

    for (i, tracer) in enumerate(tracers)
        du[i] = model(Val(tracer), 0, 0, 0, t, u..., PAR)
    end

    return nothing
end

# make sure initial values are passed in right order (Z,P,N,D)
init_conditions = (P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)
tracers = required_biogeochemical_tracers(model)
u0 = [eval(:(init_conditions.$t)) for t in tracers]

# get model parameters
p = [getfield(model, f) for f in fieldnames(typeof(model))]

tspan = (0.0, 3years)

prob = ODEProblem(model_ODEs, u0, tspan, p)

sol = solve(prob, Tsit5())

# ==================================================
# Plotting
# ==================================================

# plot(sol)

# tracer order is Z,P,N,D
p = plot(sol.t, [sol.u[i][2] for i in range(1, length(sol))]; label="P")
plot!(p, sol.t, [sol.u[i][1] for i in range(1, length(sol))]; label="Z")
plot!(p, sol.t, [sol.u[i][4] for i in range(1, length(sol))]; label="D")
savefig(p, "NPZD_box_differential_equations.png")
