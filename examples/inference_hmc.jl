# =================================================================================
# This example shows how to do HMC inference with Turing.jl for the NPZD box model.
# =================================================================================

using Agate
using Agate.Library.Light

using OrdinaryDiffEq
using LinearAlgebra
using Plots
using StatsPlots
using Turing

using Oceananigans.Units
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields

const year = years = 365day

# ==================================================
# Set up NPZD model & solver timing (save-at, stop)
# ==================================================

Δt = dt = 7days
stop_time = 1year

include(joinpath("NPZD", "tracers.jl"))

# ==================================================
# Set up DifferentialEquations
# ==================================================

init_conditions = (Z=0.05, P=0.01, N=7.0, D=0.0)
# specify which parameters to do inference for (if not all)
# should we be passing values here too ?!  the user might
# not always want to stick to the BGC defaults
params = (:μ₀, :kₙ, :lᵖᵈ, :α)
tspan = (0.0, stop_time)

prob = bgc_to_ode(NPZD, cyclical_PAR(; z=-10), init_conditions, tspan, params)

# ==================================================
# Generate noisy data
# ==================================================

sol = solve(prob, Tsit5(); saveat=Δt)
data = Array(sol) + 0.025 * randn(size(Array(sol)))

# ==================================================
# Inference
# ==================================================

@model function fit_NPZD(data, prob)
    # prior distributions centered at true value and with upper limit at ~ +2sigma
    μ₀ ~ truncated(Normal(0.6989 / day, 0.1 / day); lower=0, upper=0.9 / day)
    kₙ ~ truncated(Normal(2.3868, 0.5); lower=0, upper=3.5)
    lᵖᵈ ~ truncated(Normal(0.0101 / day, 0.01 / day); lower=0, upper=0.03 / day)
    α ~ truncated(Normal(0.1953 / day, 0.05 / day); lower=0, upper=0.3 / day)

    # observation noise
    σ ~ truncated(Cauchy(0, 1); lower=0)

    # simulate model
    p = [μ₀, kₙ, lᵖᵈ, α]
    predicted = solve(prob, Tsit5(); p=p, saveat=Δt)

    # observations
    for i in 1:length(predicted)
        data[:, i] ~ MvNormal(predicted.u[i], σ^2 * I)
    end
end

model = fit_NPZD(data, prob)
chain = sample(model, NUTS(), MCMCSerial(), 1000, 2; progress=false)

# ==================================================
# Plotting
# ==================================================

traceplot = plot(chain)
savefig(traceplot, "NPZD_trace.png")

# solve the ODE for 300 randomly picked posterior samples in the chain and plot against data
dataplot = plot(; legend=false)
posterior_samples = sample(chain[[:μ₀, :kₙ, :lᵖᵈ, :α]], 300; replace=false)
for p in eachrow(Array(posterior_samples))
    sol_p = solve(prob, Tsit5(); p=p, saveat=Δt)
    # Z and P
    plot!(
        dataplot, [sol_p.u[i][1] for i in range(1, length(sol))]; alpha=0.1, color="#BBBBBB"
    )
    plot!(
        dataplot, [sol_p.u[i][2] for i in range(1, length(sol))]; alpha=0.1, color="#BBBBBB"
    )
end

# add original simulation and generated observations to plot (Z and P)
plot!(dataplot, [sol.u[i][1] for i in range(1, length(sol))]; linewidth=1, label="Z")
plot!(dataplot, [sol.u[i][2] for i in range(1, length(sol))]; linewidth=1, label="P")
scatter!(dataplot, data[1, :]; color=1, alpha=0.3, label="")
scatter!(dataplot, data[2, :]; color=2, alpha=0.3, label="")

savefig(dataplot, "NPZD_data_fit.png")
