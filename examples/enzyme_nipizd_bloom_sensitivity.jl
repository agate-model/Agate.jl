# # [Multi-parameter bloom sensitivity with Enzyme] (@id enzyme_nipizd_bloom_sensitivity_example)

# This example differentiates a bloom-scale diagnostic with respect to several
# active NiPiZD parameters at once. The active parameter vector includes a scalar
# parameter, plankton-axis vector parameters, and predator-by-prey interaction
# matrix entries. DifferentiationInterface calls Enzyme as the automatic
# differentiation backend to compute the gradient at a reference parameter vector;
# no finite perturbation is needed to define the sensitivity ranking.

using Agate
using Agate.Library.Light: CyclicalPAR
using ADTypes: AutoEnzyme
using CairoMakie
import DifferentiationInterface
using Enzyme
using LinearAlgebra: norm
using OrdinaryDiffEq: Tsit5, solve
using Printf
using SciMLBase: remake, successful_retcode
using SciMLSensitivity

using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

const NiPiZD = Agate.Models.NiPiZD
nothing #hide

# ## Static model and active parameter map

# The active parameter vector selects growth, remineralization, grazing, and
# assimilation entries that have direct ecological interpretations in the simple
# NiPiZD food web.

const BGC = NiPiZD.construct()
const TRACER_NAMES = Tuple(required_biogeochemical_tracers(BGC))
const PLANKTON_GROUPS = Agate.Introspection.plankton_groups(BGC)

function tracer_index(tracer::Symbol)
    index = findfirst(==(tracer), TRACER_NAMES)
    isnothing(index) && error("Tracer $tracer not found in $TRACER_NAMES")
    return index
end

const PHYTOPLANKTON_INDICES = Tuple(tracer_index.(PLANKTON_GROUPS.P))

const ACTIVE = Agate.Runtime.active_parameters(BGC;
    maximum_growth_rate = (:P1, :P2),
    detritus_remineralization = true,
    palatability_matrix = ((:Z1, :P1), (:Z1, :P2), (:Z2, :P1)),
    assimilation_matrix = ((:Z1, :P1),),
)
const PARAMETER_LABELS = ACTIVE.labels

const θ_REFERENCE = copy(ACTIVE.values)

const SAVEAT = collect(range(0.0, 30day; length = 31))
const TSPAN = (first(SAVEAT), last(SAVEAT))
const PAR = CyclicalPAR(-10)
nothing #hide

# These are the quick-start box-model initial concentrations. They are expanded
# into the constructed BGC tracer order.
const INITIAL_CONCENTRATIONS = (; N = 7.0, D = 0.01,
                                 Z1 = 0.01, Z2 = 0.01,
                                 P1 = 0.01, P2 = 0.1)

function initial_conditions(::Type{T}) where {T}
    return T[getproperty(INITIAL_CONCENTRATIONS, tracer) for tracer in TRACER_NAMES]
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

# ## Bloom-scale objective

# The objective is total phytoplankton biomass at the end of a short bloom
# window. Its gradient answers a local ecological question: near the reference
# parameterization, which selected rates and interaction strengths most affect
# bloom biomass?

const SENSEALG = SciMLSensitivity.GaussAdjoint(
    autojacvec = SciMLSensitivity.EnzymeVJP(),
)

function solve_primal(theta; sensealg = nothing, saveat = SAVEAT,
                      save_start = true, save_end = true, save_everystep = false)
    problem = remake(PROBLEM; p = theta)
    common_kwargs = (; abstol = 1e-8, reltol = 1e-8, verbose = false,
                     save_start, save_end, save_everystep)
    kwargs = isnothing(saveat) ? common_kwargs : (; common_kwargs..., saveat)

    return isnothing(sensealg) ? solve(problem, Tsit5(); kwargs...) :
                                 solve(problem, Tsit5(); kwargs..., sensealg)
end

function checked_solution(theta; sensealg = nothing, saveat = SAVEAT,
                          save_start = true, save_end = true, save_everystep = false)
    sol = solve_primal(theta; sensealg, saveat, save_start, save_end, save_everystep)
    successful_retcode(sol) || error("ODE solve failed with retcode $(sol.retcode)")
    return sol
end

function solve_values(theta)
    sol = checked_solution(theta; saveat = SAVEAT)
    length(sol.t) == length(SAVEAT) || error("ODE solve returned $(length(sol.t)) saved times; expected $(length(SAVEAT))")
    return Array(sol)
end

function solve_final_values(theta; sensealg = nothing)
    sol = checked_solution(theta; sensealg, saveat = nothing,
                           save_start = false, save_end = true, save_everystep = false)
    return sol.u[end]
end

function final_total_phytoplankton_values(values)
    final_index = last(axes(values, 2))
    total = zero(eltype(values))

    for i in PHYTOPLANKTON_INDICES
        total += values[i, final_index]
    end

    return total
end

function final_total_phytoplankton(theta; sensealg = nothing)
    final_values = solve_final_values(theta; sensealg)
    total = zero(eltype(final_values))

    for i in PHYTOPLANKTON_INDICES
        total += final_values[i]
    end

    return total
end

final_total_phytoplankton_adjoint(theta) = final_total_phytoplankton(theta; sensealg = SENSEALG)

# ## Enzyme gradient through DifferentiationInterface

# DifferentiationInterface provides a common gradient API for several Julia AD
# backends. This example uses `AutoEnzyme`, so Enzyme is still the reverse-mode
# backend, while the Agate objective remains an ordinary `theta -> scalar`
# function.
#
# The scaled ranking `|θᵢ ∂J/∂θᵢ|` puts parameters with different units on a more
# comparable local-relative scale. It is a local sensitivity at `θ_REFERENCE`,
# not a finite-perturbation experiment. The differentiated solve only saves the
# final state; the denser trajectory below is for plotting the reference run.

const AD_BACKEND = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse))

enzyme_gradient(theta) = DifferentiationInterface.gradient(final_total_phytoplankton_adjoint,
                                                           AD_BACKEND,
                                                           theta)

function print_parameter_table(objective, theta, gradient, scaled_sensitivities, order)
    total = sum(scaled_sensitivities)

    @printf("Bloom sensitivity at θ_REFERENCE:\n")
    @printf("  objective final(total phytoplankton): %.8e\n", objective)
    @printf("  gradient norm:          %.8e\n", norm(gradient))
    @printf("\n  %-34s %14s %14s %14s\n", "parameter", "θ", "∂J/∂θ", "|θ ∂J/∂θ|")

    for n in order
        share = iszero(total) ? 0.0 : 100 * scaled_sensitivities[n] / total
        @printf("  %-34s %.8e % .8e %.8e  (%5.1f%%)\n",
                PARAMETER_LABELS[n], theta[n], gradient[n], scaled_sensitivities[n], share)
    end
end

function tracer_series(values, tracer::Symbol)
    return values[tracer_index(tracer), :]
end

function summed_series(values, tracers)
    total = copy(tracer_series(values, first(tracers)))

    for tracer in Iterators.drop(tracers, 1)
        total .+= tracer_series(values, tracer)
    end

    return total
end

function aggregated_trajectories(values)
    return (
        N = tracer_series(values, :N),
        D = tracer_series(values, :D),
        Z = summed_series(values, PLANKTON_GROUPS.Z),
        P = summed_series(values, PLANKTON_GROUPS.P),
    )
end

function save_diagnostic_plots(reference_values, scaled_sensitivities, order)
    days = SAVEAT ./ day
    output_directory = @__DIR__

    trajectory_path = joinpath(output_directory, "enzyme_nipizd_bloom_trajectory.png")
    sensitivity_path = joinpath(output_directory, "enzyme_nipizd_bloom_scaled_sensitivities.png")

    trajectories = aggregated_trajectories(reference_values)
    panels = (("Nutrient", trajectories.N),
              ("Detritus", trajectories.D),
              ("Total zooplankton", trajectories.Z),
              ("Total phytoplankton", trajectories.P))

    trajectory_fig = Figure(size = (1000, 720), fontsize = 14)
    for (n, (label, trajectory)) in enumerate(panels)
        row = fld(n - 1, 2) + 1
        col = mod(n - 1, 2) + 1
        ax = Axis(trajectory_fig[row, col],
                  ylabel = label,
                  xlabel = row == 2 ? "Time (days)" : "",
                  title = label)
        lines!(ax, days, trajectory, linewidth = 3)
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
              title = "Scaled bloom sensitivities",
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

reference_values = solve_values(θ_REFERENCE)
reference_objective = final_total_phytoplankton_values(reference_values)
gradient = enzyme_gradient(copy(θ_REFERENCE))
scaled_sensitivities = abs.(θ_REFERENCE .* gradient)
sensitivity_order = sortperm(scaled_sensitivities; rev = true)

print_parameter_table(reference_objective, θ_REFERENCE, gradient, scaled_sensitivities, sensitivity_order)
save_diagnostic_plots(reference_values, scaled_sensitivities, sensitivity_order)
