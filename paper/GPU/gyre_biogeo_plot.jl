using Oceananigans
using CairoMakie

function finite_colorrange(A)
    vals = vec(A)
    finite = vals[isfinite.(vals)]

    isempty(finite) && return nothing

    lo, hi = extrema(finite)

    if lo == hi
        ϵ = max(abs(lo), one(lo)) * 1e-6
        lo -= ϵ
        hi += ϵ
    end

    return (lo, hi)
end

function plot_bgc_last_timestep(filename::AbstractString;
                                tracers = (:P1, :P2, :Z1, :Z2, :N, :D),
                                slice::Symbol = :surface,
                                k::Union{Nothing, Int} = nothing,
                                j::Union{Nothing, Int} = nothing,
                                savepath::Union{Nothing, AbstractString} = nothing)

    basename = endswith(filename, ".jld2") ? filename[1:end-5] : filename
    ts = Dict(tr => FieldTimeSeries(basename, String(tr)) for tr in tracers)

    nt = length(first(values(ts)).times)
    t_last_days = first(values(ts)).times[nt] / 86400

    sample_field = first(values(ts))[nt]
    xnodes, ynodes, znodes = nodes(sample_field)

    x = collect(xnodes) ./ 1e3
    y = collect(ynodes) ./ 1e3
    z = collect(znodes)

    if slice == :surface
        k = isnothing(k) ? length(z) : k
        zlabel = "z = $(round(z[k]; digits=1)) m"
    elseif slice == :xz
        j = isnothing(j) ? cld(length(y), 2) : j
        ylabel = "y = $(round(y[j]; digits=1)) km"
    else
        error("slice must be :surface or :xz")
    end

    fig = Figure(size = (1600, 900), fontsize = 16)

    for (n, tr) in enumerate(tracers)
        row = fld1(n, 3)
        col = mod1(n, 3)

        ax = Axis(fig[row, 2col - 1])

        field = ts[tr][nt]

        if slice == :surface
            data = Array(interior(field, :, :, k))
            cr = finite_colorrange(data)

            cr === nothing && error("Tracer $(tr) has no finite values on the final :surface slice. The simulation likely produced NaNs.")

            hm = heatmap!(ax, x, y, data';
                          colorrange = cr,
                          nan_color = :transparent)

            ax.xlabel = "x (km)"
            ax.ylabel = "y (km)"
            ax.title = "$(tr) | t = $(round(t_last_days; digits=1)) d | $zlabel"
            ax.aspect = DataAspect()
        else
            data = Array(interior(field, :, j, :))
            cr = finite_colorrange(data)

            cr === nothing && error("Tracer $(tr) has no finite values on the final :xz slice. The simulation likely produced NaNs.")

            hm = heatmap!(ax, x, z, data';
                          colorrange = cr,
                          nan_color = :transparent)

            ax.xlabel = "x (km)"
            ax.ylabel = "z (m)"
            ax.title = "$(tr) | t = $(round(t_last_days; digits=1)) d | $ylabel"
        end

        Colorbar(fig[row, 2col], hm, label = string(tr))
    end

    if !isnothing(savepath)
        save(savepath, fig)
    end

    return fig
end

fig = plot_bgc_last_timestep("double_gyre_NiPiZD_noCC_noFW";
                             slice = :surface,
                             savepath = "bgc_last_surface.png")
display(fig)