# # [Reverse-mode AD sensitivity] (@id reverse_mode_ad_onyx_sensitivity_example)

# This example differentiates the final total phytoplankton biomass in the
# default Onyx model with respect to several active parameters at once. The
# active parameter vector includes a scalar parameter, plankton-axis vector
# parameters, and predator-by-prey interaction matrix entries. We then
# differentiate the same objective with respect to two grazing traits that
# jointly determine the palatability matrix.

using Agate
using Onyx
using Agate.Construction: active_traits, realize_trait_parameters
using Agate.Library.Light: CyclicalPAR
using ADTypes: AutoEnzyme
using CairoMakie
import DifferentiationInterface
using Enzyme
using OrdinaryDiffEq: Tsit5, solve
using SciMLBase: remake, successful_retcode
using SciMLSensitivity

using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Oceananigans.Units: day

nothing #hide

# ## Static model and active parameter map

# To run a model sensitivity analysis, we first choose the parameters that we
# want to vary: the active parameters.
#
# Here we vary parameters controlling phytoplankton growth, detritus recycling,
# grazing preferences, and grazing efficiency. We then measure how changes in
# each parameter affect the final total phytoplankton biomass.

const RECIPE = last(Onyx.construct_plus_recipe(; scalar_type = Float64))
const TRAITS = active_traits(Onyx.OnyxFamily(), RECIPE;
    scalar_type = Float64,
    trait_selection = (
        specificity = (:Z_1,),
        optimum_predator_prey_ratio = (:Z_1,),
    ),
)
const BGC = TRAITS.bgc
const TRACER_NAMES = Tuple(required_biogeochemical_tracers(BGC))
const PLANKTON_PFTS = Agate.Introspection.pfts(BGC)

function tracer_index(tracer::Symbol)
    index = findfirst(==(tracer), TRACER_NAMES)
    isnothing(index) && error("Tracer $tracer not found in $TRACER_NAMES")
    return index
end

const PHYTOPLANKTON = PLANKTON_PFTS.P
const PHYTOPLANKTON_INDICES = tracer_index.(PHYTOPLANKTON)

const ACTIVE = Agate.Runtime.active_parameters(BGC;
    maximum_growth_rate = (:P_1, :P_2),
    detritus_remineralization = true,
    palatability_matrix = ((:Z_1, :P_1), (:Z_1, :P_2), (:Z_2, :P_1), (:Z_2, :P_2)),
    assimilation_matrix = ((:Z_1, :P_1),),
)
const PARAMETER_LABELS = ACTIVE.labels

const θ_REFERENCE = copy(ACTIVE.values)

const SAVEAT = collect(range(0.0, 30day; length = 31))
const TSPAN = (first(SAVEAT), last(SAVEAT))
const PAR = CyclicalPAR(-10)
nothing #hide

const INITIAL_CONCENTRATIONS = (; N = 7.0, D = 0.01,
                                 Z_1 = 0.01, Z_2 = 0.01,
                                 P_1 = 0.01, P_2 = 0.1)

function initial_conditions(::Type{T}) where {T}
    return T[getproperty(INITIAL_CONCENTRATIONS, tracer) for tracer in TRACER_NAMES]
end

function onyx_problem(theta)
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

const PROBLEM = onyx_problem(θ_REFERENCE)
nothing #hide

# ## Endpoint objective

# We measure sensitivity using the final total phytoplankton biomass after a
# short model run. The gradient tells us how strongly each active parameter
# affects that final biomass near the reference parameter values.

const SENSEALG = SciMLSensitivity.GaussAdjoint(
    autojacvec = SciMLSensitivity.EnzymeVJP(),
)

function solve_problem_values(problem; sensealg = nothing, saveat = SAVEAT,
                              save_start = true, save_end = true, save_everystep = false,
                              tol = 1e-8)
    kwargs = (; abstol = tol, reltol = tol,
              saveat, save_start, save_end, save_everystep)

    if isnothing(sensealg)
        sol = solve(problem, Tsit5(); kwargs...)
        @assert successful_retcode(sol) "ODE solve failed: $(sol.retcode)"
        @assert last(sol.t) == TSPAN[end] "ODE did not reach the requested endpoint"
    else
        sol = solve(problem, Tsit5(); kwargs..., sensealg)
    end

    return Array(sol)
end

solve_values(theta; kwargs...) = solve_problem_values(remake(PROBLEM; p = theta); kwargs...)

function solve_final_values(theta; sensealg = nothing)
    values = solve_values(theta;
                          sensealg,
                          saveat = TSPAN[end],
                          save_start = false,
                          save_end = true,
                          save_everystep = false)

    return values[:, end]
end

function total_phytoplankton(values)
    total = zero(eltype(values))

    for i in PHYTOPLANKTON_INDICES
        total += values[i]
    end

    return total
end

function final_total_phytoplankton(theta; sensealg = nothing)
    return total_phytoplankton(solve_final_values(theta; sensealg))
end

final_total_phytoplankton_adjoint(theta) = final_total_phytoplankton(theta; sensealg = SENSEALG)

nothing #hide

# ## Enzyme gradient through DifferentiationInterface

# DifferentiationInterface provides a common gradient API for several Julia AD
# backends. This example uses `AutoEnzyme`, which uses Enzyme.jl as the reverse-mode
# backend.
#
# Parameters can have different units and magnitudes, so raw gradients are not
# always easy to compare. We rank sensitivities visually by `|θᵢ ∂J/∂θᵢ|`,
# which estimates how much each parameter contributes relative to its reference
# value.

const AD_BACKEND = AutoEnzyme(; mode = Enzyme.set_runtime_activity(Enzyme.Reverse))

enzyme_gradient(theta) = DifferentiationInterface.gradient(final_total_phytoplankton_adjoint,
                                                           AD_BACKEND,
                                                           theta)

function diagnostic_plots(reference_values, scaled_sensitivities, order)
    days = SAVEAT ./ day

    trajectory_fig = Figure()
    ax = Axis(trajectory_fig[1, 1],
              xlabel = "Time (days)",
              ylabel = "Total phytoplankton")
    phytoplankton = vec(sum(reference_values[PHYTOPLANKTON_INDICES, :]; dims = 1))
    lines!(ax, days, phytoplankton)

    plot_order = reverse(order)

    sensitivity_fig = Figure(size = (800, 400))
    ax = Axis(sensitivity_fig[1, 1],
              xlabel = "|θᵢ ∂J/∂θᵢ|",
              yticks = (1:length(plot_order), collect(PARAMETER_LABELS)[plot_order]))
    barplot!(ax, scaled_sensitivities[plot_order]; direction = :x)

    return (; trajectory = trajectory_fig, scaled_sensitivities = sensitivity_fig)
end

reference_values = solve_values(θ_REFERENCE)
gradient = enzyme_gradient(copy(θ_REFERENCE))
scaled_sensitivities = abs.(θ_REFERENCE .* gradient)
sensitivity_order = sortperm(scaled_sensitivities; rev = true)

plots = diagnostic_plots(reference_values, scaled_sensitivities, sensitivity_order)

# The reference trajectory shows the simulated total phytoplankton biomass used
# for the endpoint sensitivity calculation.

plots.trajectory

# The scaled sensitivities rank active parameters by `|θᵢ ∂J/∂θᵢ|`.

plots.scaled_sensitivities

# ## Construction-trait sensitivity

# The predator's optimum predator-to-prey diameter ratio determines the
# location of its prey-size preference peak, while specificity controls how
# sharply preference declines away from that peak. Both traits influence
# the palatability matrix. The default ratio of 10 is a scientific model
# assumption that we can now perturb and differentiate.
#
# `active_traits` identifies every applicable runtime entry influenced by
# the selected traits. We include those entries in the runtime parameter
# selection above, alongside the other independent parameters. This lets
# both sensitivity calculations use the same compiled ODE adjoint. Match
# the entries by label to preserve the trait map's ordering.

const τ_REFERENCE = copy(TRAITS.values)
const TRAIT_PARAMETER_INDICES = [
    only(findall(==(label), PARAMETER_LABELS)) for label in TRAITS.active_parameters.labels
]

function trait_runtime_parameters(tau)
    p = realize_trait_parameters(TRAITS, tau)
    theta = copy(θ_REFERENCE)
    theta[TRAIT_PARAMETER_INDICES] = p
    return theta
end

@assert trait_runtime_parameters(τ_REFERENCE) ≈ θ_REFERENCE

# For a trait vector τ, the runtime parameters are p = f(τ). We apply the
# chain rule in two stages: first obtain the gradient of the endpoint
# objective with respect to the runtime parameters, then propagate those
# weights through the numerical trait map using Enzyme. The weights are
# recomputed at each trait vector and held constant during the pullback.

function weighted_trait_parameters(tau, parameter_map, weights)
    p = parameter_map(tau)
    total = zero(eltype(tau))
    for i in eachindex(p, weights)
        total += p[i] * weights[i]
    end
    return total
end

function trait_gradient(tau)
    theta = trait_runtime_parameters(tau)
    @assert isfinite(final_total_phytoplankton(theta))
    weights = enzyme_gradient(theta)[TRAIT_PARAMETER_INDICES]
    @assert all(isfinite, weights)
    g = zeros(eltype(tau), length(tau))
    Enzyme.autodiff(
        Enzyme.set_runtime_activity(Enzyme.Reverse),
        Enzyme.Const(weighted_trait_parameters), Enzyme.Active,
        Enzyme.Duplicated(copy(tau), g), Enzyme.Const(TRAITS.parameter_map),
        Enzyme.Const(weights),
    )
    return g
end

# ### Verify the complete gradient

# Independently reconstruct Onyx at each perturbed trait vector. Check
# that the resulting palatability matrix agrees entry by entry with the
# runtime mapping, then solve the ordinary model with the same initial
# conditions and light forcing. This verifies the entire trait-to-objective
# gradient, including dependencies realized during model construction.

function ordinary_trait_problem(tau)
    ordinary = Onyx.construct(; scalar_type = Float64, parameters = (
        specificity = (Z_1 = tau[1],),
        optimum_predator_prey_ratio = (Z_1 = tau[2],),
    ))
    @assert Tuple(required_biogeochemical_tracers(ordinary)) == TRACER_NAMES
    mapped = Agate.Runtime.parameterized(BGC, trait_runtime_parameters(tau);
        active_parameters = ACTIVE)
    matrix = ordinary.parameters.palatability_matrix
    @assert all(isapprox(matrix[i, j], mapped.parameters.palatability_matrix[i, j];
                         rtol = 1e-12, atol = 1e-14)
                for i in axes(matrix, 1), j in axes(matrix, 2))
    @assert all(name -> name === :palatability_matrix ||
        isequal(getproperty(ordinary.parameters, name), getproperty(BGC.parameters, name)),
        keys(BGC.parameters))
    return Agate.Runtime.ode_problem(ordinary, initial_conditions(Float64), TSPAN;
        auxiliary = (; PAR = t -> PAR(Float64(t))))
end

function ordinary_trait_endpoint(tau)
    values = solve_problem_values(ordinary_trait_problem(tau);
        saveat = TSPAN[end], save_start = false, tol = 1e-10)
    endpoint = total_phytoplankton(values[:, end])
    @assert isapprox(endpoint, final_total_phytoplankton(trait_runtime_parameters(tau));
                     rtol = 1e-6, atol = 1e-8)
    return endpoint
end

function finite_difference_trait_gradient(tau)
    return map(eachindex(tau)) do i
        h = 1e-3 * max(abs(tau[i]), 0.1)
        plus, minus = copy(tau), copy(tau)
        plus[i] += h
        minus[i] -= h
        (ordinary_trait_endpoint(plus) - ordinary_trait_endpoint(minus)) / (2h)
    end
end

@assert isapprox(reference_values, solve_problem_values(ordinary_trait_problem(τ_REFERENCE));
                 rtol = 1e-6, atol = 1e-8)

# Verify the signed gradients at the default traits and at a second
# trait vector, then repeat the first gradient to check that the
# numerical mapping and derivative weights have no stale state.

τ_changed = 1.1 .* τ_REFERENCE
trait_gradients = map((τ_REFERENCE, τ_changed)) do tau
    g = trait_gradient(tau)
    fd = finite_difference_trait_gradient(tau)
    @assert all(isfinite, g) "Nonfinite trait gradient at $tau: $g"
    @assert all(isapprox(a, b; rtol = 1e-2, atol = 1e-6) for (a, b) in zip(g, fd)) "Gradient mismatch at $tau: AD = $g, finite difference = $fd"
    (; traits = copy(tau), gradient = g, finite_difference = fd, scaled = tau .* g)
end

reference_trait_gradient = first(trait_gradients).gradient
@assert all(abs.(reference_trait_gradient) .> 1e-5)
@assert isapprox(trait_gradient(τ_REFERENCE), reference_trait_gradient; rtol = 1e-8, atol = 1e-10) "Reference gradient changed after updating traits"
@assert trait_runtime_parameters(τ_REFERENCE) ≈ θ_REFERENCE

for result in trait_gradients
    println("Trait values: ", result.traits)
    for (label, g, fd, scaled) in zip(TRAITS.labels, result.gradient,
                                      result.finite_difference, result.scaled)
        println(label, ": AD = ", g, ", finite difference = ", fd,
                ", scaled sensitivity = ", scaled)
    end
end

# ### Compare trait sensitivities

# Signed scaled sensitivities, τᵢ ∂J/∂τᵢ, show how a small proportional
# increase in each trait changes the endpoint biomass near the reference
# configuration. The sign indicates whether that change increases or
# decreases the biomass. These local sensitivities can guide subsequent
# trait calibration or exploration of alternative prey-size preferences.

trait_fig = Figure(size = (800, 250))
ax = Axis(trait_fig[1, 1],
          xlabel = "τᵢ ∂J/∂τᵢ",
          yticks = (1:length(TRAITS.labels), collect(TRAITS.labels)))
barplot!(ax, first(trait_gradients).scaled; direction = :x)
trait_fig
