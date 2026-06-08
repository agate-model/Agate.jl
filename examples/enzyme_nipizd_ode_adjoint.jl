# # [Reverse-mode ODE sensitivity with Enzyme] (@id enzyme_nipizd_ode_adjoint_example)

# This example computes the gradient of a scalar objective with respect to a NiPiZD
# parameter using a reverse-mode ODE sensitivity. The ODE problem is built with
# Agate.jl's `ode_problem` helper so the biological model is constructed once and
# the active parameter enters through the ODE `p` vector.

using Agate
using Enzyme
using OrdinaryDiffEq: Tsit5, solve
using SciMLBase: remake
using SciMLSensitivity
using CairoMakie
using Printf

using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

const NiPiZD = Agate.Models.NiPiZD
const TRACERS = (:N, :D, :Z1, :Z2, :P1, :P2)
nothing #hide

# ## Static model and active parameter map

# The BGC object is static. The active value is supplied through the ODE parameter
# vector `p`, where `p[1]` is the maximum growth rate of `P1`.

const BGC = NiPiZD.construct()

required_biogeochemical_tracers(BGC) == TRACERS ||
    error("Unexpected NiPiZD tracer order: $(required_biogeochemical_tracers(BGC))")

const ACTIVE = Agate.Runtime.active_parameters(BGC; maximum_growth_rate = (:P1,))
const SAVEAT = collect(range(0.0, 365day; length=366))
const TSPAN = (first(SAVEAT), last(SAVEAT))
nothing #hide

function initial_conditions(::Type{T}) where {T}
    return T[7.0, 1.0, 0.05, 0.05, 0.01, 0.01]
end

constant_PAR(::Type{T}) where {T} = T(100.0)

function nipizd_problem(theta)
    T = eltype(theta)

    return Agate.Runtime.ode_problem(
        BGC,
        initial_conditions(T),
        TSPAN;
        p = theta,
        active_parameters = ACTIVE,
        auxiliary = (; PAR = constant_PAR(T)),
    )
end

const μ₀ = 0.7 / day
const θ₀ = [μ₀]
const PROBLEM = nipizd_problem(θ₀)
nothing #hide

# ## Scalar objective

# Reverse-mode sensitivities are most useful for scalar objectives. Here the
# objective is the final `P1` concentration after one year.

const SENSEALG = SciMLSensitivity.GaussAdjoint(
    autojacvec = SciMLSensitivity.EnzymeVJP(),
)

function p1_final_objective(theta; sensealg = nothing)
    problem = remake(PROBLEM; p = theta)

    if isnothing(sensealg)
        sol = solve(problem, Tsit5();
                    saveat = SAVEAT,
                    abstol = 1e-10,
                    reltol = 1e-10)
    else
        sol = solve(problem, Tsit5();
                    saveat = SAVEAT,
                    abstol = 1e-10,
                    reltol = 1e-10,
                    sensealg)
    end

    return sol[5, end]
end

p1_final_objective_adjoint(theta) = p1_final_objective(theta; sensealg = SENSEALG)

baseline_objective = p1_final_objective(θ₀)

# ## Enzyme gradient and finite-difference check

# We call Enzyme directly on the scalar objective. The objective closure and its
# captured static model are marked constant; only `theta` is active. The finite-
# difference check uses the primal solve without `sensealg`, so it does not run
# adjoint machinery just to estimate the reference slope.

function enzyme_gradient(theta)
    gradient = zeros(eltype(theta), length(theta))

    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
                    Enzyme.Const(p1_final_objective_adjoint),
                    Enzyme.Active,
                    Enzyme.Duplicated(theta, gradient))

    return gradient
end

function finite_difference_gradient(theta, delta)
    plus = copy(theta)
    minus = copy(theta)
    plus[1] += delta
    minus[1] -= delta

    objective_plus = p1_final_objective(plus)
    objective_minus = p1_final_objective(minus)

    return [(objective_plus - objective_minus) / (2delta)], objective_plus, objective_minus
end

function finite_difference_scan(theta; relative_steps = (1e-3, 3e-4, 1e-4, 3e-5, 1e-5))
    results = map(relative_steps) do relative_step
        delta = abs(theta[1]) * relative_step
        gradient, objective_plus, objective_minus = finite_difference_gradient(theta, delta)
        return (; relative_step, delta, gradient = gradient[1], objective_plus, objective_minus)
    end

    selected = first(filter(result -> isfinite(result.gradient) &&
                                      isfinite(result.objective_plus) &&
                                      isfinite(result.objective_minus),
                            results))

    return selected, results
end

finite_difference_sensitivity, finite_difference_results = finite_difference_scan(θ₀)
enzyme_sensitivity = enzyme_gradient(copy(θ₀))

@printf("Final P1 objective: %.8e\n", baseline_objective)
@printf("d final P1 / dμ₁, Enzyme adjoint:     %.8e\n", enzyme_sensitivity[1])
@printf("d final P1 / dμ₁, finite difference: %.8e  (relative step %.1e)\n",
        finite_difference_sensitivity.gradient,
        finite_difference_sensitivity.relative_step)
@printf("Absolute sensitivity difference: %.8e\n",
        abs(enzyme_sensitivity[1] - finite_difference_sensitivity.gradient))

# ## Plotting the trajectory

# For context, we also plot the baseline `P1` trajectory used by the objective.

sol = solve(PROBLEM, Tsit5(); saveat=SAVEAT, abstol=1e-10, reltol=1e-10)
time_days = SAVEAT ./ day
p1 = Array(sol)[5, :]

fig = Figure(; size=(900, 420), fontsize=14)

ax = Axis(fig[1, 1];
          xlabel = "time (days)",
          ylabel = "P1 concentration",
          title = "NiPiZD P1 trajectory for reverse-mode objective")
lines!(ax, time_days, p1; label="P1", linewidth=3)
axislegend(ax; position=:rb)

output_path = joinpath(@__DIR__, "enzyme_nipizd_ode_adjoint.png")
save(output_path, fig)

fig
