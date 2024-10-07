# This example shows how to integrate a Biogeochemistry (BGC) box model with DifferentialEquations.

using Agate
using DifferentialEquations
using Plots

using Oceananigans.Units
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields

const year = years = 365day

# ==================================================
# PAR
# ==================================================

const z = -10
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
PAR_f(t) = PAR⁰(t) * exp(0.2z)

# ==================================================
# Define BGC model (NPZD)
# ==================================================

# NPZD `parameters` and `tracers` are defined here
include("../examples/NPZD/model_definition.jl")

BGC = create_bgc_struct(:BGC, parameters)
add_bgc_methods(
    BGC,
    tracers,
    # we are only expecting PAR as an auxiliary field
    auxiliary_fields=[:PAR,],
    helper_functions="NPZD/functions.jl"
    )
model = BGC()

# ==================================================
# DifferentialEquations
# ==================================================

function model_DE(du, u, p, t)

    model = BGC(p...)

    PAR = PAR_f(t)

    for (i, tracer) in enumerate(tracers)
        du[i] = model(Val(tracer),0,0,0,t,u...,PAR)
    end

    return nothing
end

# make sure initial values are passed in right order (Z,P,N,D)
init_conditions = (N = 7.0, P = 0.01, Z = 0.05, D=0.0)
tracers = required_biogeochemical_tracers(model)
u0 = [eval(:(init_conditions.$t)) for t in tracers]

p = [getfield(model,f) for f in fieldnames(typeof(model))]

tspan = (0.0, 3years)

prob = ODEProblem(model_DE, u0, tspan, p)

sol = solve(prob, Tsit5())

# ==================================================
# Plotting
# ==================================================

# plot(sol)

# tracer order is Z,P,N,D
p = plot(sol.t, [sol.u[i][2] for i in range(1, length(sol))], label="P")
plot!(p, sol.t, [sol.u[i][1] for i in range(1, length(sol))], label="Z")
plot!(p, sol.t, [sol.u[i][4] for i in range(1, length(sol))], label="D")
savefig(p, "NPZD_box_differential_equations.png")
