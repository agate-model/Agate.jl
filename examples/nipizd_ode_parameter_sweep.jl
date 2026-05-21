# # [ODE parameter sweep] (@id nipizd_ode_parameter_sweep_example)

# This example shows how to use [`Agate.Diagnostics.ode_problem`](@ref) as a
# small 0D backend for parameter sweeps. The same `ODEProblem` can be solved
# directly, remade with a different constructed Agate model, or wrapped in a
# SciML `EnsembleProblem` so that each trajectory runs an independent parameter
# combination.
#
# Here we sweep over the maximum equivalent spherical diameter used to build the
# phytoplankton and zooplankton size spectra. Each grid point constructs a new
# NiPiZD biogeochemistry object, solves the same initial-value problem, and
# records the final total phytoplankton biomass.

# ## Loading dependencies

using Agate
using Agate.Diagnostics: ode_problem
using Agate.Introspection: tracer_names
using CairoMakie
using Oceananigans.Units
using OrdinaryDiffEq: Tsit5, solve
using SciMLBase: EnsembleProblem, EnsembleThreads, remake

nothing #hide

# ## Build one ODE problem

# The helper below keeps the tracer count fixed while changing the size range
# used to construct each plankton community. Keeping the tracer count fixed is a
# useful default for parameter sweeps because every trajectory has the same state
# vector layout.

function build_bgc(; phyto_max_esd=10.0, zoo_max_esd=100.0)
    phyto_size_structure = (n=2, min_esd=1.0, max_esd=phyto_max_esd, splitting=:log_splitting)
    zoo_size_structure = (n=2, min_esd=10.0, max_esd=zoo_max_esd, splitting=:log_splitting)

    return Agate.Models.NiPiZD.construct(; phyto_size_structure, zoo_size_structure)
end

# The initial conditions are passed as a `NamedTuple`. The ODE wrapper converts
# them to a state vector using Agate's tracer order, which is returned below as
# metadata for plotting and diagnostics.

initial = (N=7.0, D=0.01, Z1=0.01, Z2=0.01, P1=0.01, P2=0.01)
base_bgc = build_bgc()
metadata = ode_problem(base_bgc; initial, tspan=(0.0, 365days), return_metadata=true)
base_problem = metadata.problem
tracers = metadata.tracers

println("tracer order = ", tracers)

nothing #hide

# We can solve the base case directly.

base_solution = solve(base_problem, Tsit5(); saveat=0:30days:365days,
                      abstol=1e-6, reltol=1e-6)

println("base final state = ", base_solution.u[end])

nothing #hide

# ## Sweep a grid of parameter values

# Each ensemble trajectory receives an integer index. We map that index onto a
# pair of size-spectrum settings, build a new Agate model, and use
# `remake(problem; p=bgc_i)` to swap that model into the ODE problem. This is the
# same pattern used in calibration or Bayesian inference workflows, including
# Turing.jl models.

phyto_max_esd_grid = [6.0, 8.0, 10.0, 14.0, 18.0]
zoo_max_esd_grid = [60.0, 80.0, 100.0, 140.0, 180.0]
parameter_grid = collect(Iterators.product(phyto_max_esd_grid, zoo_max_esd_grid))

function prob_func(problem, i, repeat)
    phyto_max_esd, zoo_max_esd = parameter_grid[i]
    bgc_i = build_bgc(; phyto_max_esd, zoo_max_esd)
    return remake(problem; p=bgc_i)
end

ensemble_problem = EnsembleProblem(base_problem; prob_func)
ensemble_solution = solve(ensemble_problem, Tsit5(), EnsembleThreads();
                          trajectories=length(parameter_grid),
                          saveat=0:30days:365days,
                          abstol=1e-6,
                          reltol=1e-6)

nothing #hide

# ## Summarize the ensemble

# The output states use the same tracer order as the base problem. For this
# parameter-space plot, we extract the final total phytoplankton biomass,
# `P1 + P2`, at each grid point.

phyto_indices = findall(tracer -> startswith(string(tracer), "P"), tracers)

function final_total_phyto(solution)
    final_state = solution.u[end]
    return sum(final_state[i] for i in phyto_indices)
end

final_phyto = [final_total_phyto(ensemble_solution[i]) for i in eachindex(parameter_grid)]
final_phyto_matrix = reshape(final_phyto, length(phyto_max_esd_grid), length(zoo_max_esd_grid))

best_index = argmax(final_phyto)
best_phyto_max_esd, best_zoo_max_esd = parameter_grid[best_index]
println("maximum final total phyto = ", final_phyto[best_index])
println("at phyto max ESD = ", best_phyto_max_esd, ", zoo max ESD = ", best_zoo_max_esd)

nothing #hide

# ## Plot the parameter space

# The heatmap is a compact phase-space-style view of the sweep: one axis varies
# the phytoplankton size spectrum, the other varies the zooplankton size
# spectrum, and the color reports the final total phytoplankton biomass.

fig = Figure(; size=(780, 560), fontsize=16)
ax = Axis(fig[1, 1];
          xlabel="phytoplankton max ESD",
          ylabel="zooplankton max ESD",
          title="NiPiZD ODE parameter sweep")

hm = heatmap!(ax, phyto_max_esd_grid, zoo_max_esd_grid, final_phyto_matrix)
scatter!(ax, [best_phyto_max_esd], [best_zoo_max_esd]; marker=:star5, markersize=18)
Colorbar(fig[1, 2], hm; label="final P1 + P2")

save("nipizd_ode_parameter_sweep.png", fig; px_per_unit=1)

fig
