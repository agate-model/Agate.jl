# Benchmark the quick-start 0D NiPiZD box model through OceanBioME/Oceananigans
# against two direct ODE solves of the same Agate RHS:
#   1. a hand-written Oceananigans-style quasi-AB2 loop, and
#   2. an out-of-the-box fixed-step SSPRK22 solve from OrdinaryDiffEq.
#
# Tsit5 is timed both as a high-accuracy reference and as a faster relaxed-tolerance solve.
#
# The benchmark also writes a residual plot for P1 relative to Tsit5 to:
#   benchmark_P1_residual_vs_Tsit5.png

using Agate
using Agate.Library.Light
using Agate.Introspection: tracer_names
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using OrdinaryDiffEq: ODEProblem, Tsit5, SSPRK22, solve
using CairoMakie
using Printf
using Statistics

const NiPiZD = Agate.Models.NiPiZD
const N_PHYTO = 10
const N_ZOO = 10
const INITIAL_N = 7.0
const INITIAL_D = 0.01
const INITIAL_P = 0.01
const INITIAL_Z = 0.01
const Δt = 240minutes
const stop_time = 1095days
const QUASI_AB2_CHI = 0.1
const TRAJECTORY_TRACER = :P1
const RESIDUAL_PLOT = "benchmark_P1_residual_vs_Tsit5.png"
const TSIT5_REFERENCE_ABSTOL = 1e-10
const TSIT5_REFERENCE_RELTOL = 1e-10
const TSIT5_FAST_ABSTOL = 1e-6
const TSIT5_FAST_RELTOL = 1e-6

function timed_repeats(f, repeats)
    times = Float64[]
    allocations = Int[]
    bytes = Int[]

    for _ in 1:repeats
        GC.gc()
        result = @timed f()
        push!(times, result.time)
        push!(allocations, result.gcstats.allocd)
        push!(bytes, result.bytes)
    end

    return (; min=minimum(times), median=median(times), mean=mean(times), allocations, bytes)
end

function build_bgc(; n_phyto=N_PHYTO, n_zoo=N_ZOO)
    phyto_size_structure = (n=n_phyto, min_esd=1, max_esd=10, splitting=:log_splitting)
    zoo_size_structure = (n=n_zoo, min_esd=10, max_esd=100, splitting=:log_splitting)
    return NiPiZD.construct(; phyto_size_structure, zoo_size_structure)
end

ordered_tracers(bgc) = Tuple(tracer_names(bgc))

function initial_value(tracer)
    name = string(tracer)

    if tracer === :N
        return INITIAL_N
    elseif tracer === :D
        return INITIAL_D
    elseif startswith(name, "P")
        return INITIAL_P
    elseif startswith(name, "Z")
        return INITIAL_Z
    else
        throw(ArgumentError("No initial value rule for tracer $tracer."))
    end
end

initial_conditions(tracers) = (; (tracer => initial_value(tracer) for tracer in tracers)...)

function build_box_model(; n_phyto=N_PHYTO, n_zoo=N_ZOO)
    grid = BoxModelGrid()
    light_attenuation = FunctionFieldPAR(; grid)
    bgc = build_bgc(; n_phyto, n_zoo)
    tracers = ordered_tracers(bgc)
    bgc_model = Biogeochemistry(bgc; light_attenuation)
    model = BoxModel(; biogeochemistry=bgc_model)
    set!(model; initial_conditions(tracers)...)
    return (; model, tracers)
end

read_box_state(model, tracers) = [getproperty(model.fields, tracer).data[1, 1, 1] for tracer in tracers]

function run_oceanbiome_box(; n_phyto=N_PHYTO, n_zoo=N_ZOO, collect_daily=false)
    (; model, tracers) = build_box_model(; n_phyto, n_zoo)
    nsteps = Int(round(stop_time / Δt))

    if collect_daily
        samples = Vector{Vector{Float64}}()
        sample_interval = Int(round(day / Δt))
        push!(samples, read_box_state(model, tracers))
        for step in 1:nsteps
            time_step!(model, Δt)
            if step % sample_interval == 0
                push!(samples, read_box_state(model, tracers))
            end
        end
        return samples
    else
        for _ in 1:nsteps
            time_step!(model, Δt)
        end
        return read_box_state(model, tracers)
    end
end

function build_direct_rhs(; n_phyto=N_PHYTO, n_zoo=N_ZOO)
    bgc = build_bgc(; n_phyto, n_zoo)
    tracers = ordered_tracers(bgc)

    function rhs!(du, u, t)
        PAR = CyclicalPAR(-10)(t)
        for (i, tracer) in enumerate(tracers)
            du[i] = bgc(Val(tracer), 0, 0, 0, t, u..., PAR)
        end
        return nothing
    end

    u0 = [initial_value(tracer) for tracer in tracers]
    return (; rhs!, u0, tracers)
end

function build_ode_problem(; n_phyto=N_PHYTO, n_zoo=N_ZOO)
    (; rhs!, u0, tracers) = build_direct_rhs(; n_phyto, n_zoo)

    function ode_rhs!(du, u, _, t)
        rhs!(du, u, t)
        return nothing
    end

    return (; problem=ODEProblem(ode_rhs!, u0, (0.0, stop_time)), tracers)
end

function run_tsit5(; n_phyto=N_PHYTO, n_zoo=N_ZOO, collect_daily=false,
                    abstol=TSIT5_REFERENCE_ABSTOL,
                    reltol=TSIT5_REFERENCE_RELTOL)
    (; problem) = build_ode_problem(; n_phyto, n_zoo)

    if collect_daily
        solution = solve(problem, Tsit5(); saveat=0:day:stop_time,
                         abstol, reltol)
        return [Vector(u) for u in solution.u]
    else
        solution = solve(problem, Tsit5(); save_everystep=false,
                         abstol, reltol)
        return Vector(solution.u[end])
    end
end

run_tsit5_reference(; kwargs...) = run_tsit5(; kwargs...,
                                             abstol=TSIT5_REFERENCE_ABSTOL,
                                             reltol=TSIT5_REFERENCE_RELTOL)

run_tsit5_fast(; kwargs...) = run_tsit5(; kwargs...,
                                        abstol=TSIT5_FAST_ABSTOL,
                                        reltol=TSIT5_FAST_RELTOL)

function run_direct_quasi_ab2(; n_phyto=N_PHYTO, n_zoo=N_ZOO, collect_daily=false, chi=QUASI_AB2_CHI)
    (; rhs!, u0, tracers) = build_direct_rhs(; n_phyto, n_zoo)
    nsteps = Int(round(stop_time / Δt))
    u = copy(u0)
    f_previous = similar(u)
    f_current = similar(u)
    samples = collect_daily ? Vector{Vector{Float64}}() : nothing
    sample_interval = Int(round(day / Δt))

    rhs!(f_previous, u, 0.0)

    if collect_daily
        push!(samples, copy(u))
    end

    for step in 1:nsteps
        if step == 1
            @. u = u + Δt * f_previous
        else
            @. u = u + Δt * ((1.5 + chi) * f_current - (0.5 + chi) * f_previous)
            f_previous, f_current = f_current, f_previous
        end

        rhs!(f_current, u, step * Δt)

        if collect_daily && step % sample_interval == 0
            push!(samples, copy(u))
        end
    end

    return collect_daily ? samples : u
end

function run_ssprk22(; n_phyto=N_PHYTO, n_zoo=N_ZOO, collect_daily=false)
    (; problem) = build_ode_problem(; n_phyto, n_zoo)

    if collect_daily
        solution = solve(problem, SSPRK22(); dt=Δt, adaptive=false,
                         saveat=0:day:stop_time, save_everystep=false)
        return [Vector(u) for u in solution.u]
    else
        solution = solve(problem, SSPRK22(); dt=Δt, adaptive=false,
                         save_everystep=false)
        return Vector(solution.u[end])
    end
end

function max_daily_difference(samples_a, samples_b)
    n = min(length(samples_a), length(samples_b))
    return maximum(maximum(abs.(samples_a[i] .- samples_b[i])) for i in 1:n)
end

trajectory(samples, tracer_index) = [sample[tracer_index] for sample in samples]

function plot_tracer_residuals(box_samples, quasi_ab2_samples, ssprk22_samples, tsit5_fast_samples, tsit5_reference_samples, tracers;
                               tracer=TRAJECTORY_TRACER,
                               filename=RESIDUAL_PLOT)
    tracer_index = findfirst(==(tracer), tracers)
    tracer_index === nothing && throw(ArgumentError("Tracer $tracer not found in $tracers."))

    n = min(length(box_samples), length(quasi_ab2_samples), length(ssprk22_samples),
            length(tsit5_fast_samples), length(tsit5_reference_samples))
    days = collect(0:(n - 1))
    box_values = trajectory(box_samples[1:n], tracer_index)
    quasi_ab2_values = trajectory(quasi_ab2_samples[1:n], tracer_index)
    ssprk22_values = trajectory(ssprk22_samples[1:n], tracer_index)
    tsit5_fast_values = trajectory(tsit5_fast_samples[1:n], tracer_index)
    tsit5_reference_values = trajectory(tsit5_reference_samples[1:n], tracer_index)

    box_residual = box_values .- tsit5_reference_values
    quasi_ab2_residual = quasi_ab2_values .- tsit5_reference_values
    ssprk22_residual = ssprk22_values .- tsit5_reference_values
    tsit5_fast_residual = tsit5_fast_values .- tsit5_reference_values

    fig = Figure(size=(900, 520))
    ax = Axis(fig[1, 1],
              xlabel="time (days)",
              ylabel="$(tracer) residual",
              title="$(tracer) residuals relative to Tsit5")

    lines!(ax, days, box_residual, label="OceanBioME BoxModel - Tsit5")
    lines!(ax, days, quasi_ab2_residual, label="Direct quasi-AB2 χ=$(QUASI_AB2_CHI) - Tsit5")
    lines!(ax, days, ssprk22_residual, label="DiffEq SSPRK22 - Tsit5 reference")
    lines!(ax, days, tsit5_fast_residual, label="Tsit5 fast - Tsit5 reference")
    hlines!(ax, [0.0], linestyle=:dash, label="Tsit5 reference baseline")
    axislegend(ax, position=:rt)
    save(filename, fig)
    return filename
end

function print_stats(label, stats)
    @printf("%-28s min %.4f s | median %.4f s | mean %.4f s | bytes/run %.3e\n",
            label, stats.min, stats.median, stats.mean, median(stats.bytes))
end

function main(; repeats=5, n_phyto=N_PHYTO, n_zoo=N_ZOO)
    bgc = build_bgc(; n_phyto, n_zoo)
    tracers = ordered_tracers(bgc)
    initial = initial_conditions(tracers)

    println("Quick-start NiPiZD 0D benchmark")
    println("n_phyto = $n_phyto, n_zoo = $n_zoo, n_tracers = $(length(tracers))")
    println("tracer order = $tracers")
    println("initial conditions = $initial")
    println("Δt = $(Δt / minute) minutes, stop_time = $(stop_time / day) days, repeats = $repeats")
    println("direct quasi-AB2 χ = $QUASI_AB2_CHI")
    println("Tsit5 reference tolerances = (abstol=$(TSIT5_REFERENCE_ABSTOL), reltol=$(TSIT5_REFERENCE_RELTOL))")
    println("Tsit5 fast tolerances = (abstol=$(TSIT5_FAST_ABSTOL), reltol=$(TSIT5_FAST_RELTOL))")
    println()

    run_oceanbiome_box(; n_phyto, n_zoo)
    run_direct_quasi_ab2(; n_phyto, n_zoo)
    run_ssprk22(; n_phyto, n_zoo)
    run_tsit5_fast(; n_phyto, n_zoo)
    run_tsit5_reference(; n_phyto, n_zoo, collect_daily=true)

    box_stats = timed_repeats(() -> run_oceanbiome_box(; n_phyto, n_zoo), repeats)
    quasi_ab2_stats = timed_repeats(() -> run_direct_quasi_ab2(; n_phyto, n_zoo), repeats)
    ssprk22_stats = timed_repeats(() -> run_ssprk22(; n_phyto, n_zoo), repeats)
    tsit5_fast_stats = timed_repeats(() -> run_tsit5_fast(; n_phyto, n_zoo), repeats)
    tsit5_reference_stats = timed_repeats(() -> run_tsit5_reference(; n_phyto, n_zoo), repeats)

    box_samples = run_oceanbiome_box(; n_phyto, n_zoo, collect_daily=true)
    quasi_ab2_samples = run_direct_quasi_ab2(; n_phyto, n_zoo, collect_daily=true)
    ssprk22_samples = run_ssprk22(; n_phyto, n_zoo, collect_daily=true)
    tsit5_fast_samples = run_tsit5_fast(; n_phyto, n_zoo, collect_daily=true)
    tsit5_reference_samples = run_tsit5_reference(; n_phyto, n_zoo, collect_daily=true)

    box_quasi_ab2_difference = max_daily_difference(box_samples, quasi_ab2_samples)
    box_ssprk22_difference = max_daily_difference(box_samples, ssprk22_samples)
    box_tsit5_difference = max_daily_difference(box_samples, tsit5_reference_samples)
    quasi_ab2_tsit5_difference = max_daily_difference(quasi_ab2_samples, tsit5_reference_samples)
    ssprk22_tsit5_difference = max_daily_difference(ssprk22_samples, tsit5_reference_samples)
    tsit5_fast_reference_difference = max_daily_difference(tsit5_fast_samples, tsit5_reference_samples)
    plot_file = plot_tracer_residuals(box_samples, quasi_ab2_samples, ssprk22_samples,
                                      tsit5_fast_samples, tsit5_reference_samples, tracers)

    print_stats("OceanBioME BoxModel", box_stats)
    print_stats("Direct quasi-AB2 χ=0.1", quasi_ab2_stats)
    print_stats("DiffEq SSPRK22", ssprk22_stats)
    print_stats("Tsit5 fast", tsit5_fast_stats)
    print_stats("Tsit5 reference", tsit5_reference_stats)
    @printf("%-28s %.6e\n", "BoxModel vs qAB2 daily diff", box_quasi_ab2_difference)
    @printf("%-28s %.6e\n", "BoxModel vs SSPRK22 diff", box_ssprk22_difference)
    @printf("%-28s %.6e\n", "BoxModel vs Tsit5 daily diff", box_tsit5_difference)
    @printf("%-28s %.6e\n", "qAB2 vs Tsit5 daily diff", quasi_ab2_tsit5_difference)
    @printf("%-28s %.6e\n", "SSPRK22 vs Tsit5 diff", ssprk22_tsit5_difference)
    @printf("%-28s %.3fx\n", "BoxModel / qAB2 median", box_stats.median / quasi_ab2_stats.median)
    @printf("%-28s %.3fx\n", "BoxModel / SSPRK22 median", box_stats.median / ssprk22_stats.median)
    println("$(TRAJECTORY_TRACER) residual plot = $plot_file")

    return (; box_stats, quasi_ab2_stats, ssprk22_stats, tsit5_fast_stats, tsit5_reference_stats,
            box_quasi_ab2_difference, box_ssprk22_difference, box_tsit5_difference,
            quasi_ab2_tsit5_difference, ssprk22_tsit5_difference, tsit5_fast_reference_difference,
            n_phyto, n_zoo, tracers, initial, plot_file)
end

main()
