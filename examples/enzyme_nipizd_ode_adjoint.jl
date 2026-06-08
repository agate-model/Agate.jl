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

const ACTIVE_PARAMETERS = (; maximum_growth_rate = (; P1 = 1))
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
        active_parameters = ACTIVE_PARAMETERS,
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

function p1_final_objective(theta)
    problem = remake(PROBLEM; p = theta)
    sol = solve(problem, Tsit5();
                saveat = SAVEAT,
                abstol = 1e-10,
                reltol = 1e-10,
                sensealg = SENSEALG)

    return sol[5, end]
end

baseline_objective = p1_final_objective(θ₀)

# ## Enzyme gradient and finite-difference check

# We call Enzyme directly on the scalar objective. The objective closure and its
# captured static model are marked constant; only `theta` is active.

function enzyme_gradient(theta)
    gradient = zeros(eltype(theta), length(theta))

    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
                    Enzyme.Const(p1_final_objective),
                    Enzyme.Active,
                    Enzyme.Duplicated(theta, gradient))

    return gradient
end

function finite_difference_gradient(theta, delta)
    plus = copy(theta)
    minus = copy(theta)
    plus[1] += delta
    minus[1] -= delta
    return [(p1_final_objective(plus) - p1_final_objective(minus)) / (2delta)]
end

δ = μ₀ * 1e-6
enzyme_sensitivity = enzyme_gradient(copy(θ₀))
finite_difference_sensitivity = finite_difference_gradient(θ₀, δ)

@printf("Final P1 objective: %.8e\n", baseline_objective)
@printf("d final P1 / dμ₁, Enzyme adjoint:     %.8e\n", enzyme_sensitivity[1])
@printf("d final P1 / dμ₁, finite difference: %.8e\n", finite_difference_sensitivity[1])
@printf("Absolute sensitivity difference: %.8e\n", abs(enzyme_sensitivity[1] - finite_difference_sensitivity[1]))

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
