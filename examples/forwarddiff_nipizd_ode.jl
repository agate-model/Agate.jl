# # [Forward-mode AD sensitivity] (@id forwarddiff_nipizd_ode_example)

# In this example we use ForwardDiff.jl to compute the sensitivity of P1 in the N2P2ZD model to the maximum growth rate parameter.

# ## Loading dependencies
# The example uses Agate.jl for the ecosystem model, OrdinaryDiffEq.jl for a small standalone ODE solve,
# ForwardDiff.jl for the derivative, and CairoMakie.jl for plotting.

using Agate
using ForwardDiff
using OrdinaryDiffEq: ODEProblem, Tsit5, solve
using CairoMakie
using Printf

using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

const NiPiZD = Agate.Models.NiPiZD
const TRACERS = (:N, :D, :Z1, :Z2, :P1, :P2)
nothing #hide

# ## Constructing an AD-active model

# Agate uses an explicit constructor-boundary scalar type contract.
# For ordinary Oceananigans.jl simulations this scalar type comes from the grid, or defaults to `Float64`.
# For ForwardDiff.jl, we pass the active scalar type explicitly with `scalar_type = typeof(mu)`.

function nipizd_model_with_p1_growth_rate(mu)
    T = typeof(mu)

    return NiPiZD.construct(;
        scalar_type=T,
        parameters=(; maximum_growth_rate=[zero(T), zero(T), mu, T(0.7 / day)]),
    )
end
nothing #hide

# Here `mu` is the maximum growth rate of `P1`.
# The two leading zeros correspond to zooplankton entries.
# The second phytoplankton growth rate is held fixed.

# ## Standalone ODE problem

# We define a small standalone ODE that calls Agate's generated biological tendency functions directly.

function initial_conditions(::Type{T}) where {T}
    return T[7.0, 1.0, 0.05, 0.05, 0.01, 0.01]
end

constant_PAR(::Type{T}) where {T} = T(100.0)

function solve_nipizd(mu; saveat=range(0.0, 365day; length=366))
    T = typeof(mu)
    bgc = nipizd_model_with_p1_growth_rate(mu)

    required_biogeochemical_tracers(bgc) == TRACERS ||
        error("Unexpected NiPiZD tracer order: $(required_biogeochemical_tracers(bgc))")

    function rhs!(du, u, _, t)
        PAR = constant_PAR(T)
        for (i, tracer) in enumerate(TRACERS)
            du[i] = bgc(Val(tracer), 0, 0, 0, t, u..., PAR)
        end
        return nothing
    end

    u0 = initial_conditions(T)
    problem = ODEProblem(rhs!, u0, (first(saveat), last(saveat)))
    return solve(problem, Tsit5(); saveat=saveat, abstol=1e-10, reltol=1e-10)
end
nothing #hide

# We expose the `P1` trajectory (e.g. biomass values over time) as a vector-valued function of one parameter.
# This is the function that ForwardDiff differentiates.

function p1_trajectory(theta; saveat=range(0.0, 365day; length=366))
    sol = solve_nipizd(theta[1]; saveat=saveat)
    values = reduce(hcat, sol.u)
    return vec(values[5, :])
end

function p1_solution(mu; saveat=range(0.0, 365day; length=366))
    sol = solve_nipizd(mu; saveat=saveat)
    values = reduce(hcat, sol.u)
    return vec(values[5, :])
end
nothing #hide

# ## ForwardDiff and finite differences

# We compute the time-dependent sensitivity of `P1` to its own maximum growth rate.
# A central finite difference provides a simple independent check.

function finite_difference_p1_trajectory(mu0, delta; saveat)
    plus = p1_trajectory([mu0 + delta]; saveat=saveat)
    minus = p1_trajectory([mu0 - delta]; saveat=saveat)
    return (plus .- minus) ./ (2delta)
end

saveat = collect(range(0.0, 365day; length=366))
mu0 = 0.7 / day
delta = mu0 * 1e-6

baseline = p1_solution(mu0; saveat=saveat)
J = ForwardDiff.jacobian(theta -> p1_trajectory(theta; saveat=saveat), [mu0])
forwarddiff_sensitivity = J[:, 1]
finite_difference_sensitivity = finite_difference_p1_trajectory(mu0, delta; saveat=saveat)

@printf("Final dP1/dmu, ForwardDiff:       %.8e\n", forwarddiff_sensitivity[end])
@printf("Final dP1/dmu, finite difference: %.8e\n", finite_difference_sensitivity[end])
@printf(
    "Maximum absolute sensitivity difference: %.8e\n",
    maximum(abs.(forwarddiff_sensitivity .- finite_difference_sensitivity)),
)

# ## Plotting

# The top panel shows P1 biomass over time.
# The bottom panel compares the ForwardDiff sensitivity with the central finite-difference estimate.

time_days = saveat ./ day
fig = Figure(; size=(900, 620), fontsize=14)

ax1 = Axis(fig[1, 1]; xlabel="time (days)", ylabel="P1 concentration", title="NiPiZD P1 biomass")
lines!(ax1, time_days, baseline; label="P1", linewidth=3)
axislegend(ax1; position=:rb)

ax2 = Axis(fig[2, 1]; xlabel="time (days)", ylabel="sensitivity to P1 growth rate", title="ForwardDiff sensitivity vs finite-difference sensitivity")
lines!(ax2, time_days, forwarddiff_sensitivity; label="dP1/dμ₁, ForwardDiff", linewidth=3)
lines!(ax2, time_days, finite_difference_sensitivity; linestyle=:dash, label="dP1/dμ₁, finite difference", linewidth=3)
axislegend(ax2; position=:rb)

output_path = joinpath(@__DIR__, "forwarddiff_nipizd_ode_sensitivity.png")
save(output_path, fig)

fig