# # [Multi-parameter reverse-mode ODE sensitivity with Enzyme] (@id enzyme_nipizd_multiparameter_adjoint_example)

# This example computes the gradient of a scalar trajectory misfit with respect
# to several NiPiZD parameters at once. The active parameter vector includes a
# scalar parameter, plankton-axis vector parameters, and predator-by-prey
# interaction matrix entries. The biological model is constructed once and the
# active values enter only through the ODE `p` vector. The setup follows the
# quick-start box-model defaults: plausible NiPiZD parameter values,
# quick-start initial tracer concentrations, and the default seasonal PAR curve
# used by `FunctionFieldPAR` in a box model.

using Agate
using Agate.Library.Light: CyclicalPAR
using CairoMakie
using Enzyme
using LinearAlgebra: norm
using OrdinaryDiffEq: Tsit5, solve
using Printf
using SciMLBase: remake, successful_retcode
using SciMLSensitivity

using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

const NiPiZD = Agate.Models.NiPiZD
const TRACERS = (:N, :D, :Z1, :Z2, :P1, :P2)
nothing #hide

# ## Static model and active parameter map

# The BGC object uses the package defaults. The active parameter vector supplies
# selected values from several parameter domains:
#
# - `maximum_growth_rate.P1` and `.P2` are phytoplankton growth-rate entries.
# - `detritus_remineralization` controls recycling from detritus back to nutrient.
# - `palatability_matrix` and `assimilation_matrix` entries control grazing
#   interactions.

const BGC = NiPiZD.construct()

required_biogeochemical_tracers(BGC) == TRACERS ||
    error("Unexpected NiPiZD tracer order: $(required_biogeochemical_tracers(BGC))")

# `ACTIVE` names the active leaves and carries the corresponding reference
# values from the default BGC.
const ACTIVE = Agate.Runtime.active_parameters(BGC;
    maximum_growth_rate = (:P1, :P2),
    detritus_remineralization = true,
    palatability_matrix = ((:Z1, :P1), (:Z1, :P2), (:Z2, :P1)),
    assimilation_matrix = ((:Z1, :P1),),
)
const PARAMETER_LABELS = ACTIVE.labels
const θ_REFERENCE = ACTIVE.values

const θ_INITIAL = [
    1.75e-5,
    9.60e-6,
    1.40e-6,
    0.95,
    0.31,
    0.12,
    0.30,
]

const SAVEAT = collect(range(0.0, 365day; length=366))
const TSPAN = (first(SAVEAT), last(SAVEAT))
const PAR = CyclicalPAR(-10)
nothing #hide

# These are the quick-start box-model initial concentrations, reordered to match
# the BGC tracer order `(:N, :D, :Z1, :Z2, :P1, :P2)`.
function initial_conditions(::Type{T}) where {T}
    return T[7.0, 0.01, 0.01, 0.01, 0.01, 0.1]
end

function nipizd_problem(theta)
    T = eltype(theta)

    return Agate.Runtime.ode_problem(
        BGC,
        initial_conditions(T),
        TSPAN;
        p = theta,
        active_parameters = ACTIVE,
        auxiliary = (; PAR = t -> PAR(T(t))),
    )
end

const PROBLEM = nipizd_problem(θ_REFERENCE)
nothing #hide

# ## Synthetic target trajectory

# The target is generated from `θ_REFERENCE`. The objective asks: if a plausible
# initial parameter vector `θ_INITIAL` differs from this reference, which active
# parameters most strongly reduce the trajectory mismatch back toward the
# reference quick-start ecosystem trajectory?

function solve_primal(theta; sensealg = nothing)
    problem = remake(PROBLEM; p = theta)
    kwargs = (; saveat = SAVEAT, abstol = 1e-10, reltol = 1e-10, verbose = false)
    return isnothing(sensealg) ? solve(problem, Tsit5(); kwargs...) :
                                 solve(problem, Tsit5(); kwargs..., sensealg)
end

function solve_values(theta; sensealg = nothing)
    sol = solve_primal(theta; sensealg)

    successful_retcode(sol) || error("ODE solve failed with retcode $(sol.retcode)")
    length(sol.t) == length(SAVEAT) || error("ODE solve returned $(length(sol.t)) saved times; expected $(length(SAVEAT))")

    return Array(sol)
end

const TARGET_SOLUTION = solve_primal(θ_REFERENCE)
successful_retcode(TARGET_SOLUTION) || error("Target solve failed with retcode $(TARGET_SOLUTION.retcode)")
const TARGET_VALUES = Array(TARGET_SOLUTION)

const θ₀ = θ_INITIAL
const θ₀_SOLUTION = solve_primal(θ₀)
successful_retcode(θ₀_SOLUTION) || error("Initial perturbed solve failed with retcode $(θ₀_SOLUTION.retcode)")
length(θ₀_SOLUTION.t) == length(SAVEAT) || error("Initial perturbed solve returned $(length(θ₀_SOLUTION.t)) saved times; expected $(length(SAVEAT))")
nothing #hide

# ## Scalar objective

# Reverse-mode sensitivities are most useful for scalar objectives. Here the
# objective is a concentration-space trajectory mismatch across all six tracers.
# Each tracer residual is scaled by its target trajectory mean, so the nutrient
# inventory does not swamp plankton or detritus differences purely because it has
# larger absolute concentration units.

const SENSEALG = SciMLSensitivity.GaussAdjoint(
    autojacvec = SciMLSensitivity.EnzymeVJP(),
)

const OBJECTIVE_TRACERS = ((:N, 1), (:D, 2), (:Z1, 3), (:Z2, 4), (:P1, 5), (:P2, 6))
const CONCENTRATION_SCALE_FLOOR = 1e-8

function tracer_scale(target, idx)
    trajectory = target[idx, :]
    return sum(abs, trajectory) / length(trajectory) + CONCENTRATION_SCALE_FLOOR
end

function scaled_concentration_misfit(values, target, idx)
    scale = tracer_scale(target, idx)
    residual = (values[idx, :] .- target[idx, :]) ./ scale
    return sum(abs2, residual) / length(residual)
end

function trajectory_misfit(theta; sensealg = nothing)
    values = solve_values(theta; sensealg)
    total = zero(eltype(values))

    for (_, idx) in OBJECTIVE_TRACERS
        total += scaled_concentration_misfit(values, TARGET_VALUES, idx)
    end

    return total / length(OBJECTIVE_TRACERS)
end

trajectory_misfit_adjoint(theta) = trajectory_misfit(theta; sensealg = SENSEALG)

# ## Enzyme gradient

# We call Enzyme directly on the scalar objective. This multiparameter example is
# intentionally Enzyme-only: finite-difference sanity checks are kept in the
# single-parameter example and package tests, where they are easier to keep
# stable and interpretable. Long nonlinear ecological integrations can be very
# sensitive to finite-difference perturbations, especially for interaction
# parameters.

function enzyme_gradient(theta)
    gradient = zeros(eltype(theta), length(theta))

    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse),
                    Enzyme.Const(trajectory_misfit_adjoint),
                    Enzyme.Active,
                    Enzyme.Duplicated(theta, gradient))

    return gradient
end

function print_parameter_table(reference, initial)
    @printf("Active parameter values:\n")
    @printf("  %-34s %14s %14s\n", "parameter", "reference", "initial θ₀")
    for n in eachindex(PARAMETER_LABELS)
        @printf("  %-34s %.8e %.8e\n", PARAMETER_LABELS[n], reference[n], initial[n])
    end
end

function print_sensitivity_summary(theta, gradient)
    scaled_sensitivities = abs.(theta .* gradient)
    total = sum(scaled_sensitivities)
    order = sortperm(scaled_sensitivities; rev = true)

    @printf("\nScaled sensitivity ranking |θᵢ ∂J/∂θᵢ|:\n")
    for n in order
        share = iszero(total) ? 0.0 : 100 * scaled_sensitivities[n] / total
        @printf("  %-34s %.8e  (%5.1f%%)\n", PARAMETER_LABELS[n], scaled_sensitivities[n], share)
    end

    return scaled_sensitivities, order
end

function save_diagnostic_plots(current_values, scaled_sensitivities, order)
    days = SAVEAT ./ day
    output_directory = @__DIR__

    trajectory_path = joinpath(output_directory, "enzyme_nipizd_multiparameter_trajectories.png")
    sensitivity_path = joinpath(output_directory, "enzyme_nipizd_multiparameter_scaled_sensitivities.png")

    trajectory_fig = Figure(size = (1100, 850), fontsize = 14)
    for (n, (label, idx)) in enumerate(OBJECTIVE_TRACERS)
        row = fld(n - 1, 2) + 1
        col = mod(n - 1, 2) + 1
        ax = Axis(trajectory_fig[row, col],
                  ylabel = string(label),
                  xlabel = row == 3 ? "Time (days)" : "",
                  title = "$(label)")
        lines!(ax, days, TARGET_VALUES[idx, :], label = "reference target", linewidth = 3)
        lines!(ax, days, current_values[idx, :], label = "initial θ₀", linewidth = 3)
        n == 1 && axislegend(ax, position = :rt)
    end
    save(trajectory_path, trajectory_fig)

    plot_order = reverse(order)
    labels = collect(PARAMETER_LABELS)[plot_order]
    scores = scaled_sensitivities[plot_order]
    positive_scores = filter(!iszero, scores)
    plot_floor = isempty(positive_scores) ? 1.0 : minimum(positive_scores) / 10
    plot_scores = max.(scores, plot_floor)

    sensitivity_fig = Figure(size = (1000, 520), fontsize = 14)
    ax = Axis(sensitivity_fig[1, 1],
              title = "Scaled parameter sensitivities",
              xlabel = "|θᵢ ∂J/∂θᵢ|",
              xscale = log10,
              yticks = (1:length(labels), labels))
    barplot!(ax, plot_scores; direction = :x)
    save(sensitivity_path, sensitivity_fig)

    @printf("\nSaved diagnostic plots:\n")
    @printf("  %s\n", trajectory_path)
    @printf("  %s\n", sensitivity_path)

    return trajectory_path, sensitivity_path
end

baseline_loss = trajectory_misfit(θ₀)
gradient = enzyme_gradient(copy(θ₀))
current_values = solve_values(θ₀)

@printf("Trajectory misfit at θ₀: %.8e\n", baseline_loss)
print_parameter_table(θ_REFERENCE, θ₀)
@printf("Gradient norm: %.8e\n", norm(gradient))
@printf("\nSelected Enzyme gradient entries:\n")
for n in eachindex(PARAMETER_LABELS)
    @printf("  %-34s % .8e\n", PARAMETER_LABELS[n], gradient[n])
end

scaled_sensitivities, sensitivity_order = print_sensitivity_summary(θ₀, gradient)
save_diagnostic_plots(current_values, scaled_sensitivities, sensitivity_order)
