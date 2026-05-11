# plot_last_year_10P10Z.jl
#
# Post-process a saved Oceananigans / Agate NiPiZD JLD2 output file.
# Writes:
#   1. annual mean over the final simulation year, surface/top-20-m x-y maps
#   2. annual mean over the final simulation year, x-mean y-z transects
#
# Edit the paths below if your output file or directory names differ.

using Oceananigans
using Oceananigans.Units: day
using JLD2
using Statistics
using Printf

ENV["GKSwstype"] = "100"
using Plots

gr()

const FT = Float64
const year = 365day

const output_file = "./paper/GPU/double_gyre_N10P10ZD_from_spunup_3yr.jld2"
const output_prefix = "double_gyre_N10P10ZD_from_spunup_3yr_last_year"
const output_dir = "plots_last_year"

const Lx = 3180e3
const Ly = 3180e3
const H = 4000.0
const Nx = 60
const Ny = 60
const Nz = 30

const P_NAMES = Symbol.("P" .* string.(1:10))
const Z_NAMES = Symbol.("Z" .* string.(1:10))
const TRACER_NAMES = vcat(P_NAMES, Z_NAMES)

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

const z_faces = couespel_z_faces(FT; H=H, Nz=Nz)
const x = collect(range(-Lx / 2, Lx / 2; length=Nx)) ./ 1e3
const y = collect(range(-Ly / 2, Ly / 2; length=Ny)) ./ 1e3
const z = 0.5 .* (z_faces[1:end-1] .+ z_faces[2:end])

function open_timeseries(filepath::String, name::Symbol)
    return FieldTimeSeries(filepath, String(name); backend=OnDisk())
end

function tracer_exists(filepath::String, name::Symbol)
    try
        fts = open_timeseries(filepath, name)
        return length(fts.times) > 0
    catch err
        @warn "Skipping missing or unreadable tracer" tracer=name exception=(err, catch_backtrace())
        return false
    end
end

function final_year_window(filepath::String, names)
    first_name = first(filter(name -> tracer_exists(filepath, name), names))
    times = Float64.(open_timeseries(filepath, first_name).times)
    isempty(times) && error("No output times found in $(filepath).")

    stop_time = maximum(times)
    start_time = stop_time - year
    return start_time, stop_time
end

function time_indices_and_weights(times::Vector{Float64}, window_start::Float64, window_end::Float64)
    indices = findall(t -> window_start < t <= window_end, times)
    isempty(indices) && error("No saved snapshots found between $(window_start / day) and $(window_end / day) days.")

    selected_times = times[indices]
    weights = zeros(Float64, length(indices))

    for n in eachindex(indices)
        left = n == firstindex(indices) ? window_start : 0.5 * (selected_times[n - 1] + selected_times[n])
        right = n == lastindex(indices) ? window_end : 0.5 * (selected_times[n] + selected_times[n + 1])
        weights[n] = max(0.0, right - left)
    end

    total_weight = sum(weights)
    total_weight <= 0 && error("Computed non-positive time weights for final-year average.")
    return indices, weights ./ total_weight
end

function top_layer_weights(z_faces; top_depth::Float64=20.0)
    weights = zeros(Float64, length(z_faces) - 1)
    lower_limit = -abs(top_depth)
    upper_limit = 0.0

    for k in eachindex(weights)
        lower = z_faces[k]
        upper = z_faces[k + 1]
        weights[k] = max(0.0, min(upper, upper_limit) - max(lower, lower_limit))
    end

    sum(weights) <= 0 && error("No grid-cell overlap found in the top $(top_depth) m.")
    return weights ./ sum(weights)
end

function annual_mean_field(filepath::String, name::Symbol, window_start::Float64, window_end::Float64)
    fts = open_timeseries(filepath, name)
    times = Float64.(fts.times)
    indices, weights = time_indices_and_weights(times, window_start, window_end)

    mean_field = nothing
    for (idx, weight) in zip(indices, weights)
        snapshot = Float64.(Array(fts[idx]))
        if mean_field === nothing
            mean_field = zeros(Float64, size(snapshot))
        end
        mean_field .+= weight .* snapshot
    end

    return mean_field
end

function top20_mean_map(field::Array{Float64, 3})
    weights = top_layer_weights(z_faces; top_depth=20.0)
    surface = zeros(Float64, size(field, 1), size(field, 2))

    for k in axes(field, 3)
        weights[k] == 0 && continue
        surface .+= weights[k] .* field[:, :, k]
    end

    return surface
end

function yz_transect(field::Array{Float64, 3})
    return dropdims(mean(field, dims=1), dims=1)
end

function tracer_title(name::Symbol)
    s = String(name)
    startswith(s, "P") && return "Phyto $(s[2:end])"
    startswith(s, "Z") && return "Zoo $(s[2:end])"
    return s
end

function plot_last_year_surface_maps(filepath::String, window_start::Float64, window_end::Float64)
    plt = plot(layout=(4, 5), size=(1900, 1300),
               plot_title="Final-year annual mean, top 20 m")

    for (panel, name) in enumerate(TRACER_NAMES)
        tracer_exists(filepath, name) || continue
        @info "Computing final-year top-20-m mean" tracer=name
        field = annual_mean_field(filepath, name, window_start, window_end)
        surface = top20_mean_map(field)

        heatmap!(plt[panel], x, y, permutedims(surface),
                 xlabel="x (km)", ylabel="y (km)",
                 title=tracer_title(name),
                 colorbar_title=String(name),
                 aspect_ratio=:equal)
    end

    mkpath(output_dir)
    outfile = joinpath(output_dir, "$(output_prefix)_surface_top20m_annual_mean.png")
    savefig(plt, outfile)
    @info "Wrote $(outfile)"
    return outfile
end

function plot_last_year_yz_transects(filepath::String, window_start::Float64, window_end::Float64)
    plt = plot(layout=(4, 5), size=(1900, 1300),
               plot_title="Final-year annual mean, x-mean y-z transect")

    for (panel, name) in enumerate(TRACER_NAMES)
        tracer_exists(filepath, name) || continue
        @info "Computing final-year y-z transect" tracer=name
        field = annual_mean_field(filepath, name, window_start, window_end)
        transect = yz_transect(field)

        heatmap!(plt[panel], y, z, permutedims(transect),
                 xlabel="y (km)", ylabel="z (m)",
                 title=tracer_title(name),
                 colorbar_title=String(name))
    end

    mkpath(output_dir)
    outfile = joinpath(output_dir, "$(output_prefix)_yz_transect_annual_mean.png")
    savefig(plt, outfile)
    @info "Wrote $(outfile)"
    return outfile
end

function main()
    isfile(output_file) || error("Could not find output file: $(output_file)")

    window_start, window_end = final_year_window(output_file, TRACER_NAMES)
    @printf("Using final-year averaging window: %.3f to %.3f days\n", window_start / day, window_end / day)

    surface_file = plot_last_year_surface_maps(output_file, window_start, window_end)
    transect_file = plot_last_year_yz_transects(output_file, window_start, window_end)

    println("Wrote plots:")
    println("  ", surface_file)
    println("  ", transect_file)
    return nothing
end

main()
