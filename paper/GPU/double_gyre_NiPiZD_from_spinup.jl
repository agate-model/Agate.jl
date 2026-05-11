# Thirty-year double gyre NiPiZD run initialized from a spun-up physical state.
# The script writes summary plots after the run:
#   - domain-mean tracer time series over the full 30-year run
#   - depth profiles after year 10 and year 30
#   - surface maps after year 10 and year 30
#   - monthly-mean surface maps for the final year ending at year 10 and year 30
#
# Notes:
#   - initializes u, v, T, and S from double_gyre_physics_terminal_state.jld2
#   - restores the saved model time so the seasonal cycle continues in phase
#   - initializes NiPiZD tracers from standard defaults
#   - monthly means are computed from the saved output snapshots, so they are
#     best if out_interval_days is not too large
#
# If you prefer to start at the beginning of the archived forcing year instead of the
# terminal physical state, change physics_state_file to:
#   "double_gyre_physics_cycle_init.jld2"

using Oceananigans
using Oceananigans.Units: minute, day
using Oceananigans.Advection: UpwindBiased
using Oceananigans.OutputWriters: JLD2Writer
using Oceananigans.OutputReaders: FieldTimeSeries, OnDisk
using Oceananigans.Grids: Center, node
using Printf
using Statistics

using CUDA
CUDA.allowscalar(false)
using JLD2

using Agate
using Agate.Library.Light
using Agate.Models: NiPiZD
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers
using Adapt

using Plots

gr()

if !CUDA.functional()
    error("CUDA is not functional. Check nvidia-smi and CUDA.jl setup.")
end

const FT   = Float32
const year = 365day
const arch = GPU()

const Lx = 3180e3
const Ly = 3180e3
const H  = 4000.0
const Nx = 60
const Ny = 60
const Nz = 30

const month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
const month_labels  = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

function couespel_z_faces(::Type{FT}; H=4000.0, Nz=30) where {FT}
    Δz1 = collect(range(FT(10.0),  FT(14.0);  length=10))
    Δz2 = collect(range(FT(20.0),  FT(56.0);  length=10))
    Δz3 = collect(range(FT(100.0), FT(500.0); length=10))

    Δz1 .*= FT(120.0) / sum(Δz1)
    Δz2 .*= FT(380.0) / sum(Δz2)
    Δz3 .*= FT(H - 500.0) / sum(Δz3)

    Δz = vcat(Δz3, Δz2, Δz1)

    z_faces = zeros(FT, Nz + 1)
    z_faces[1] = -FT(H)
    for k in 1:Nz
        z_faces[k + 1] = z_faces[k] + Δz[k]
    end
    z_faces[end] = FT(0)
    return z_faces
end

z_faces = couespel_z_faces(FT; H=H, Nz=Nz)

grid = RectilinearGrid(arch, FT;
                       size=(Nx, Ny, Nz),
                       x = (-FT(Lx / 2), FT(Lx / 2)),
                       y = (-FT(Ly / 2), FT(Ly / 2)),
                       z = z_faces,
                       topology = (Bounded, Bounded, Bounded))

coriolis = BetaPlane(; latitude=FT(30),
                     radius=FT(6.371e6),
                     rotation_rate=FT(7.292115e-5))

const α  = FT(2.0e-4)
const βS = FT(7.6e-4)

eos = LinearEquationOfState(; thermal_expansion=α,
                              haline_contraction=βS)

buoyancy = SeawaterBuoyancy(; equation_of_state=eos)

const νh    = FT(1.0e5)
const κh    = FT(3.0e3)
const νz_bg = FT(1.0e-4)
const κz_bg = FT(1.0e-5)

function build_closure()
    Az = HorizontalScalarDiffusivity(; ν=νh, κ=κh)
    Kz = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization();
                                   ν=νz_bg, κ=κz_bg)
    return (Az, Kz)
end

const ρ0         = FT(1026.0)
const τ0         = FT(0.10)
const τ_season   = FT(0.015)
const h_surface  = FT(50.0)
const τT_restore = 30day
const S0         = FT(35.0)

@inline ηy(y) = (y + FT(Ly / 2)) / FT(Ly)

@inline function τx_surface(y, t)
    seasonal = τ0 + τ_season * FT(sin(2π * t / year))
    return -seasonal * FT(cos(2π * ηy(y)))
end

@inline function T_air(y, t)
    η = ηy(y)
    Tmean   = FT(15.0) + FT(10.0) * FT(cos(π * η))
    Tseason = FT(5.0)  * FT(sin(2π * t / year))
    return Tmean + Tseason
end

@inline τx_flux(x, y, t) = τx_surface(y, t) / ρ0

@inline function PAR0(x, y, t)
    seasonal = 1 - cos((t + 15day) * 2π / year)
    gaussian = exp(-((mod(t, year) - 200day) / (50day))^2)
    pulse    = 1 / (1 + 0.2 * gaussian) + 2
    return FT(60) * FT(seasonal) * FT(pulse)
end

@inline function k_bg(y)
    η = ηy(y)
    return FT(0.04) + FT(0.04) * FT(η)
end

@inline PAR_f(x, y, z, t) = PAR0(x, y, t) * exp(k_bg(y) * z)

light_attenuation = FunctionFieldPAR(; grid, PAR_f=PAR_f)

bgc_obj = NiPiZD.construct()

bgc_cpu = try
    bgc_obj()
catch
    bgc_obj
end

bgc_tracers_any = required_biogeochemical_tracers(bgc_cpu)
bgc_tracers = bgc_tracers_any isa Symbol ? (bgc_tracers_any,) : Tuple(bgc_tracers_any)

bgc_gpu = Adapt.adapt(CuArray, bgc_cpu)
biogeochemistry = Biogeochemistry(bgc_gpu; light_attenuation)

u_bcs = FieldBoundaryConditions(
    top    = FluxBoundaryCondition(τx_flux),
    bottom = FluxBoundaryCondition(FT(0)),
    south  = FluxBoundaryCondition(FT(0)),
    north  = FluxBoundaryCondition(FT(0))
)

v_bcs = FieldBoundaryConditions(
    top    = FluxBoundaryCondition(FT(0)),
    bottom = FluxBoundaryCondition(FT(0)),
    west   = FluxBoundaryCondition(FT(0)),
    east   = FluxBoundaryCondition(FT(0))
)

T_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(FT(0)),
                                bottom=FluxBoundaryCondition(FT(0)),
                                west=FluxBoundaryCondition(FT(0)),
                                east=FluxBoundaryCondition(FT(0)),
                                south=FluxBoundaryCondition(FT(0)),
                                north=FluxBoundaryCondition(FT(0)))

S_bcs = T_bcs

@inline function T_restoring_discrete(i, j, k, grid, clock, fields)
    x, y, z = node(i, j, k, grid, Center(), Center(), Center())
    T = @inbounds fields.T[i, j, k]
    return ifelse(z > -h_surface,
                  (T_air(y, clock.time) - T) / τT_restore,
                  FT(0))
end

function build_model()
    closure      = build_closure()
    momentum_adv = UpwindBiased(order=3)
    tracer_adv   = UpwindBiased(order=1)

    T_restoring = Forcing(T_restoring_discrete; discrete_form=true)
    all_tracers = Tuple(unique([:T, :S, bgc_tracers...]))

    model = HydrostaticFreeSurfaceModel(
        grid = grid,
        free_surface = ImplicitFreeSurface(),
        coriolis = coriolis,
        buoyancy = buoyancy,
        tracers = all_tracers,
        closure = closure,
        boundary_conditions = (u=u_bcs, v=v_bcs, T=T_bcs, S=S_bcs),
        momentum_advection = momentum_adv,
        tracer_advection   = tracer_adv,
        timestepper = :SplitRungeKutta3,
        forcing = (T = T_restoring,),
        biogeochemistry = biogeochemistry
    )

    return model
end

function load_physics_state!(model; filepath::String="double_gyre_physics_terminal_state.jld2")
    if !isfile(filepath)
        error("Could not find physics state file: $(filepath)")
    end

    jldopen(filepath, "r") do file
        u = FT.(file["u"])
        v = FT.(file["v"])
        T = FT.(file["T"])
        S = FT.(file["S"])

        set!(model,
             u = u,
             v = v,
             T = T,
             S = S)

        if haskey(file, "time")
            model.clock.time = Float64(file["time"])
        end
    end

    @info "Loaded spun-up physics from $(filepath)"
    return nothing
end

function initialize_bgc_tracers!(model)
    defaults = Dict{Symbol, FT}(
        :P1 => FT(0.01), :P2 => FT(0.01),
        :Z1 => FT(0.05), :Z2 => FT(0.05),
        :N  => FT(0.1),
        :D  => FT(1.0)
    )

    for (name, val) in defaults
        if haskey(model.tracers, name)
            set!(model.tracers[name], val)
        end
    end

    return nothing
end

to_cpu(a) = Array(interior(a))

function format_elapsed_time(seconds)
    total_seconds = round(Int, seconds)
    hours = total_seconds ÷ 3600
    minutes = (total_seconds % 3600) ÷ 60
    secs = total_seconds % 60
    return "$(hours) h $(minutes) min $(secs) s"
end

function tracer_label(name::Symbol)
    labels = Dict(
        :T  => "Temperature",
        :S  => "Salinity",
        :N  => "Nutrient",
        :P1 => "Phytoplankton 1",
        :P2 => "Phytoplankton 2",
        :Z1 => "Zooplankton 1",
        :Z2 => "Zooplankton 2",
        :D  => "Detritus"
    )
    return get(labels, name, String(name))
end

function month_index_from_time(t)
    day_of_year = mod(Float64(t / day), 365.0)
    cumulative = 0.0
    for (m, len) in enumerate(month_lengths)
        next_cumulative = cumulative + len
        if cumulative <= day_of_year < next_cumulative
            return m
        end
        cumulative = next_cumulative
    end
    return 12
end

function open_timeseries(filepath::String, name::Symbol)
    return FieldTimeSeries(filepath, String(name); backend=OnDisk())
end

function tracer_exists(filepath::String, name::Symbol)
    try
        fts = open_timeseries(filepath, name)
        return length(fts.times) > 0
    catch
        return false
    end
end

function snapshot_at_or_before(filepath::String, name::Symbol, target_time)
    fts = open_timeseries(filepath, name)
    times = Float64.(fts.times)
    idx = findlast(t -> t <= target_time + 1e-6, times)
    idx === nothing && error("No snapshot found for $(name) at or before target time $(target_time / day) days")
    return Array(fts[idx]), times[idx]
end

function plot_domain_mean_timeseries(filepath::String, filename::String; tracer_names=[:P1, :P2, :Z1, :Z2])
    plt = plot(xlabel="Time (years)", ylabel="Domain mean", title="Domain-mean NiPiZD tracers")

    for name in tracer_names
        tracer_exists(filepath, name) || continue
        fts = open_timeseries(filepath, name)
        n = length(fts.times)
        tyears = Float64.(fts.times) ./ year
        means = zeros(Float64, n)

        for i in 1:n
            means[i] = mean(Array(fts[i]))
        end

        plot!(plt, tyears, means, label=String(name), lw=2)
    end

    savefig(plt, "$(filename)_timeseries.png")
    return nothing
end

function plot_profile_at_time(filepath::String, target_time, filename::String, label::String; tracer_names=[:T, :N, :P1, :P2, :Z1, :Z2])
    zc = Array(znodes(grid, Center(), Center(), Center()))
    zplot = reverse(zc)

    plt = plot(xlabel="Value", ylabel="Depth (m)", title="$(label) horizontally averaged profiles")

    for name in tracer_names
        tracer_exists(filepath, name) || continue
        data, _ = snapshot_at_or_before(filepath, name, target_time)
        prof = vec(mean(data, dims=(1, 2)))
        plot!(plt, reverse(prof), zplot, label=String(name), lw=2)
    end

    savefig(plt, "$(filename)_profiles.png")
    return nothing
end

function plot_surface_map_at_time(filepath::String, name::Symbol, target_time, filename::String, label::String)
    tracer_exists(filepath, name) || return nothing

    x = Array(xnodes(grid, Center(), Center(), Center()))
    y = Array(ynodes(grid, Center(), Center(), Center()))
    field, _ = snapshot_at_or_before(filepath, name, target_time)
    data = field[:, :, end]

    plt = heatmap(x ./ 1e3, y ./ 1e3, permutedims(data);
                  xlabel="x (km)",
                  ylabel="y (km)",
                  title="$(label) surface $(tracer_label(name))",
                  colorbar_title=String(name),
                  aspect_ratio=:equal)

    savefig(plt, "$(filename)_surface_$(String(name)).png")
    return nothing
end

function compute_window_monthly_surface_means(filepath::String, name::Symbol, window_start, window_end)
    tracer_exists(filepath, name) || error("Tracer $(name) not found in $(filepath)")

    fts = open_timeseries(filepath, name)
    times = Float64.(fts.times)

    monthly_sums = [zeros(Float64, Nx, Ny) for _ in 1:12]
    monthly_counts = zeros(Int, 12)

    for i in eachindex(times)
        t = times[i]
        if window_start < t <= window_end
            m = month_index_from_time(t)
            snapshot = Array(fts[i])
            surface = snapshot[:, :, end]
            monthly_sums[m] .+= surface
            monthly_counts[m] += 1
        end
    end

    monthly_means = Vector{Matrix{Float64}}(undef, 12)
    for m in 1:12
        if monthly_counts[m] > 0
            monthly_means[m] = monthly_sums[m] ./ monthly_counts[m]
        else
            monthly_means[m] = fill(NaN, Nx, Ny)
        end
    end

    return monthly_means, monthly_counts
end

function robust_clims(monthly_means)
    finite_values = Float64[]
    for field in monthly_means
        vals = field[isfinite.(field)]
        append!(finite_values, vec(vals))
    end

    isempty(finite_values) && error("No finite values found for monthly means.")

    lo, hi = quantile(finite_values, [0.02, 0.98])
    if !isfinite(lo) || !isfinite(hi) || hi <= lo
        lo = minimum(finite_values)
        hi = maximum(finite_values)
    end
    if hi <= lo
        hi = lo + eps(Float64)
    end

    return lo, hi
end

function plot_window_monthly_surface_maps(filepath::String, name::Symbol, window_start, window_end, filename::String, label::String)
    tracer_exists(filepath, name) || return nothing

    monthly_means, monthly_counts = compute_window_monthly_surface_means(filepath, name, window_start, window_end)
    x = Array(xnodes(grid, Center(), Center(), Center())) ./ 1e3
    y = Array(ynodes(grid, Center(), Center(), Center())) ./ 1e3
    cmin, cmax = robust_clims(monthly_means)

    plt = plot(layout=(3, 4), size=(1400, 1000),
               plot_title="$(label) monthly mean surface $(tracer_label(name))")

    for m in 1:12
        subtitle = "$(month_labels[m]) (n=$(monthly_counts[m]))"
        heatmap!(plt[m], x, y, permutedims(monthly_means[m]);
                 xlabel="x (km)",
                 ylabel="y (km)",
                 title=subtitle,
                 clims=(cmin, cmax),
                 aspect_ratio=:equal)
    end

    savefig(plt, "$(filename)_monthly_surface_$(String(name)).png")
    return nothing
end

function write_plots(output_file::String, filename::String, start_time)
    isfile(output_file) || return nothing

    plot_domain_mean_timeseries(output_file, filename)

    milestones = ((10, "year10"), (30, "year30"))

    for (yr, tag) in milestones
        target_time = start_time + yr * year
        window_start = target_time - year
        prefix = "$(filename)_$(tag)"
        label = "After $(yr) years"

        plot_profile_at_time(output_file, target_time, prefix, label)

        for name in (:N, :P1, :P2, :Z1, :Z2)
            plot_surface_map_at_time(output_file, name, target_time, prefix, label)
        end

        for name in (:P1, :P2, :Z1, :Z2)
            plot_window_monthly_surface_maps(output_file, name, window_start, target_time, prefix, label)
        end
    end

    @info "Wrote milestone plots with prefix $(filename)_year10_* and $(filename)_year30_*"
    return nothing
end

function write_runtime_summary(filename::String, elapsed_seconds, simulated_days, start_time, end_time)
    simulated_years = simulated_days / 365
    speed_days_per_second = simulated_days / elapsed_seconds
    speed_years_per_day = simulated_years / (elapsed_seconds / 86400)

    summary = """
Run time summary
================
Wall-clock elapsed time : $(format_elapsed_time(elapsed_seconds))
Elapsed seconds         : $(round(elapsed_seconds, digits=2))
Simulated days          : $(round(simulated_days, digits=3))
Simulated years         : $(round(simulated_years, digits=3))
Model start time (days) : $(round(start_time / day, digits=3))
Model end time (days)   : $(round(end_time / day, digits=3))
Speed (sim days / sec)  : $(round(speed_days_per_second, digits=5))
Speed (sim years / day) : $(round(speed_years_per_day, digits=3))
"""

    open("$(filename)_runtime_summary.txt", "w") do io
        write(io, summary)
    end

    println(summary)
    @info "Wrote runtime summary to $(filename)_runtime_summary.txt"
    return nothing
end

function progress(sim, start_time)
    c = sim.model.clock
    @printf("iter: %d | elapsed: %.4f yr | absolute: %.4f yr | Δt: %.2f min\n",
            c.iteration,
            (c.time - start_time) / year,
            c.time / year,
            sim.Δt / minute)
    return nothing
end

function run_simulation(; stop_years::Float64 = 30.0,
                          Δt0 = 20minute,
                          out_interval_days::Int = 15,
                          physics_state_file::String = "double_gyre_physics_terminal_state.jld2",
                          filename::String = "double_gyre_NiPiZD_from_spunup_30yr",
                          make_plots::Bool = true)

    model = build_model()
    load_physics_state!(model; filepath=physics_state_file)
    initialize_bgc_tracers!(model)

    start_time = model.clock.time
    stop_time  = start_time + stop_years * year

    simulation = Simulation(model; Δt=Δt0, stop_time=stop_time)

    wizard = TimeStepWizard(cfl=0.5,
                            max_change=1.1,
                            min_change=0.5,
                            max_Δt=30minute)

    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
    simulation.callbacks[:progress] = Callback(sim -> progress(sim, start_time), IterationInterval(200))

    output_file = "$(filename).jld2"

    if out_interval_days > 0
        simulation.output_writers[:tracers] = JLD2Writer(
            model,
            model.tracers;
            filename = output_file,
            schedule = TimeInterval(out_interval_days * day),
            overwrite_existing = true
        )
        @info "Writing JLD2 outputs to $(output_file) every $(out_interval_days) days"
    else
        @info "Output disabled (out_interval_days <= 0)."
    end

    wall_start = time()
    run!(simulation)
    elapsed_seconds = time() - wall_start

    simulated_days = (model.clock.time - start_time) / day
    write_runtime_summary(filename, elapsed_seconds, simulated_days, start_time, model.clock.time)

    if make_plots
        write_plots(output_file, filename, start_time)
    end

    return nothing
end

run_simulation(; stop_years = 30.0,
               Δt0 = 20minute,
               out_interval_days = 15,
               physics_state_file = "double_gyre_physics_terminal_state.jld2",
               filename = "double_gyre_NiPiZD_from_spunup_30yr",
               make_plots = true)
