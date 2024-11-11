# This example shows how to integrate a Biogeochemistry (BGC) box model with DifferentialEquations.
using Agate
using Agate.Library.Light

using OrdinaryDiffEq
using Plots

using Oceananigans.Units

const year = years = 365day

# ==================================================
# Define BGC model (NPZD)
# ==================================================

include(joinpath("NPZD", "tracers.jl"))

# ==================================================
# DifferentialEquations
# ==================================================

init_conditions = (N=7.0, P=0.01, Z=0.05, D=0.0)
tspan = (0.0, 3years)

prob = bgc_to_ode(NPZD, cyclical_PAR(; z=-10), init_conditions, tspan)

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
