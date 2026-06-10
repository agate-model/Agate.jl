# # [Parameter sensitivity with reverse-mode AD] (@id ad_nipizd_parameter_sensitivity_example)

# This example differentiates the final total phytoplankton biomass in the
# default NiPiZD model with respect to several active parameters at once. The
# active parameter vector includes a scalar parameter, plankton-axis vector
# parameters, and predator-by-prey interaction matrix entries.

using Agate
using Agate.Library.Light: CyclicalPAR
using ADTypes: AutoEnzyme
using CairoMakie
import DifferentiationInterface
using Enzyme
using LinearAlgebra: norm
using OrdinaryDiffEq: Tsit5, solve
using Printf
using SciMLBase: remake
using SciMLSensitivity

using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

const NiPiZD = Agate.Models.NiPiZD
nothing #hide

# ## Static model and active parameter map

# To run a model sensitivity analysis, we first choose the parameters that we
# want to vary: the active parameters.
#
# Here we vary parameters controlling phytoplankton growth, detritus recycling,
# grazing preferences, and grazing efficiency. We then measure how changes in
# each parameter affect the final total phytoplankton biomass.

const BGC = NiPiZD.construct()
const TRACER_NAMES = Tuple(required_biogeochemical_tracers(BGC))
const PLANKTON_GROUPS = Agate.Introspection.plankton_groups(BGC)

function tracer_index(tracer::Symbol)
    index = findfirst(==(tracer), TRACER_NAMES)
    isnothing(index) && error("Tracer $tracer not found in $TRACER_NAMES")
    return index
end

const PHYTOPLANKTON = PLANKTON_GROUPS.P
const PHYTOPLANKTON_INDICES = tracer_index.(PHYTOPLANKTON)

const ACTIVE = Agate.Runtime.active_parameters(BGC;
    maximum_growth_rate = (:P1, :P2),
    detritus_remineralization = true,
    interactions = (;
        palatability = ((:Z1, :P1), (:Z1, :P2), (:Z2, :P1)),
        assimilation = ((:Z1, :P1),),
    ),
)
const PARAMETER_LABELS = ACTIVE.labels

const θ_REFERENCE = copy(ACTIVE.values)

const SAVEAT = collect(range(0.0, 30day; length = 31))
const TSPAN = (first(SAVEAT), last(SAVEAT))
const PAR = CyclicalPAR(-10)
nothing #hide

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

# ## Endpoint objective

# We measure sensitivity using the final total phytoplankton biomass after a
# short model run. The gradient tells us how strongly each active parameter
# affects that final biomass near the reference parameter values.

const SENSEALG = SciMLSensitivity.GaussAdjoint(
    autojacvec = SciMLSensitivity.EnzymeVJP(),
)

function solve_values(theta; sensealg = nothing, saveat = SAVEAT,
                      save_start = true, save_end = true, save_everystep = false)
    problem = remake(PROBLEM; p = theta)
    kwargs = (; abstol = 1e-8, reltol = 1e-8, verbose = false,
              saveat, save_start, save_end, save_everystep)

    sol = isnothing(sensealg) ? solve(problem, Tsit5(); kwargs...) :
                                solve(problem, Tsit5(); kwargs..., sensealg)

    return Array(sol)
end

function solve_final_values(theta; sensealg = nothing)
    values = solve_values(theta;
                          sensealg,
                          saveat = TSPAN[end],
                          save_start = false,
                          save_end = true,
                          save_everystep = false)

    return values[:, end]
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
# backends. This example uses `AutoEnzyme`, which uses Enzyme.jl as the reverse-mode
# backend.
#
# Parameters can have different units and magnitudes, so raw gradients are not
# always easy to compare. We rank sensitivities by `|θᵢ ∂J/∂θᵢ|`, which estimates
# how much each parameter contributes relative to its reference value.

const AD_BACKEND = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse))

enzyme_gradient(theta) = DifferentiationInterface.gradient(final_total_phytoplankton_adjoint,
                                                           AD_BACKEND,
                                                           theta)

function print_parameter_table(objective, theta, gradient, scaled_sensitivities, order)
    total = sum(scaled_sensitivities)

    @printf("Parameter sensitivity at θ_REFERENCE:\n")
    @printf("  objective final(total phytoplankton): %.8e\n", objective)
    @printf("  gradient norm:          %.8e\n", norm(gradient))
    @printf("\n  %-34s %14s %14s %14s\n", "parameter", "θ", "∂J/∂θ", "|θ ∂J/∂θ|")

    for n in order
        share = iszero(total) ? 0.0 : 100 * scaled_sensitivities[n] / total
        @printf("  %-34s %.8e % .8e %.8e  (%5.1f%%)\n",
                PARAMETER_LABELS[n], theta[n], gradient[n], scaled_sensitivities[n], share)
    end
end

function save_diagnostic_plots(reference_values, scaled_sensitivities, order)
    days = SAVEAT ./ day
    output_directory = @__DIR__

    trajectory_path = joinpath(output_directory, "ad_nipizd_parameter_sensitivity_trajectory.png")
    sensitivity_path = joinpath(output_directory, "ad_nipizd_parameter_sensitivity_scaled_sensitivities.png")

    trajectory_fig = Figure()
    ax = Axis(trajectory_fig[1, 1],
              xlabel = "Time (days)",
              ylabel = "Total phytoplankton")
    phytoplankton = vec(sum(reference_values[PHYTOPLANKTON_INDICES, :]; dims = 1))
    lines!(ax, days, phytoplankton)
    save(trajectory_path, trajectory_fig)

    plot_order = reverse(order)

    sensitivity_fig = Figure(size = (800, 400))
    ax = Axis(sensitivity_fig[1, 1],
              xlabel = "|θᵢ ∂J/∂θᵢ|",
              yticks = (1:length(plot_order), collect(PARAMETER_LABELS)[plot_order]))
    barplot!(ax, scaled_sensitivities[plot_order]; direction = :x)
    save(sensitivity_path, sensitivity_fig)

    return trajectory_path, sensitivity_path
end

reference_values = solve_values(θ_REFERENCE)
reference_objective = final_total_phytoplankton(θ_REFERENCE)
gradient = enzyme_gradient(copy(θ_REFERENCE))
scaled_sensitivities = abs.(θ_REFERENCE .* gradient)
sensitivity_order = sortperm(scaled_sensitivities; rev = true)

print_parameter_table(reference_objective, θ_REFERENCE, gradient, scaled_sensitivities, sensitivity_order)
save_diagnostic_plots(reference_values, scaled_sensitivities, sensitivity_order)
