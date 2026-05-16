# # [Predator optimal prey size] (@id predator_optimal_prey_size_example)

# This example compares the default NiPiZD palatability matrix with a model
# constructed using a different preferred predator:prey size ratio
# (`optimum_predator_prey_ratio`, abbreviated here as `Vopt`). It then runs the
# two configurations in the quickstart 0D box-model setup and compares their
# tracer trajectories.

using Agate
using Agate.Introspection: interaction_table, tracer_names
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using CairoMakie
using Printf

nothing #hide

# ## Construct the two ecosystem configurations

default_bgc = Agate.Models.NiPiZD.construct()
default_pal = interaction_table(default_bgc, :palatability)
nothing #hide

# By default, predators in the NiPiZD model prefer prey that are 10 times smaller than themselves (`Vopt` = 10). 
# Here we construct a variant with an optimum predator:prey size ratio of 5.
# NiPiZD's default plankton-tracer order is `Z1, Z2, P1, P2`. `Vopt` is a
# consumer trait, so the zooplankton entries are set to 5 and the phytoplankton
# entries remain zero.

vopt_bgc = Agate.Models.NiPiZD.construct(;
    parameters=(; optimum_predator_prey_ratio=[5.0, 5.0, 0.0, 0.0])
)
vopt_pal = interaction_table(vopt_bgc, :palatability)

nothing #hide

# ## Helper for labelled matrix plots

function plot_matrix!(fig, position, table; title)
    mat = Matrix(table.matrix)
    nrows, ncols = size(mat)

    ax = Axis(
        fig[position...];
        title,
        xlabel=string(table.column_axis),
        ylabel=string(table.row_axis),
        xticks=(1:ncols, string.(table.columns)),
        yticks=(1:nrows, string.(table.rows)),
        yreversed=true,
        aspect=DataAspect(),
    )

    hm = heatmap!(ax, 1:ncols, 1:nrows, mat'; colorrange=(0, 1))

    for row in 1:nrows, col in 1:ncols
        text!(
            ax,
            col,
            row;
            text=@sprintf("%.2f", mat[row, col]),
            align=(:center, :center),
            fontsize=12,
        )
    end

    return hm
end

nothing #hide

# ## Run both configurations in the quickstart box model

function run_box_model(bgc, filename)
    light_attenuation = FunctionFieldPAR(; grid=BoxModelGrid())
    bgc_model = Biogeochemistry(bgc; light_attenuation)
    full_model = BoxModel(; biogeochemistry=bgc_model)

    set!(full_model; N=7.0, P1=0.01, Z1=0.01, P2=0.1, Z2=0.01, D=0.01)

    simulation = Simulation(full_model; Δt=240minutes, stop_time=1095days)

    simulation.output_writers[:fields] = JLD2Writer(
        full_model,
        full_model.fields;
        filename,
        schedule=TimeInterval(1day),
        overwrite_existing=true,
    )

    run!(simulation)

    return filename
end

default_filename = run_box_model(default_bgc, "predator_optimal_prey_size_default.jld2")
vopt_filename = run_box_model(vopt_bgc, "predator_optimal_prey_size_vopt5.jld2")

nothing #hide

# ## Load the output time series

tracer_syms = tracer_names(default_bgc)

default_timeseries = (;
    (s => FieldTimeSeries(default_filename, string(s))[1, 1, 1, :] for s in tracer_syms)...
)

vopt_timeseries = (;
    (s => FieldTimeSeries(vopt_filename, string(s))[1, 1, 1, :] for s in tracer_syms)...
)

times = FieldTimeSeries(default_filename, string(first(tracer_syms))).times

nothing #hide

# ## Plot matrix and trajectory comparisons

fig = Figure(; size=(900, 1100), fontsize=16)

hm_pal_default = plot_matrix!(fig, (1, 1), default_pal; title="Vopt = 10 (default)")

plot_matrix!(fig, (1, 2), vopt_pal; title="Vopt = 5")

Colorbar(fig[1, 3], hm_pal_default; label="palatability")

Label(fig[0, 1:3], "Predator optimal prey size and box-model response"; fontsize=22)

for (idx, sym) in enumerate(tracer_syms)
    row = 2 + div(idx - 1, 2)
    col = 1 + mod(idx - 1, 2)

    ax = Axis(
        fig[row, col];
        ylabel=string(sym),
        xlabel="Days",
        title="$(sym) concentration (mmol N / m³)",
    )

    lines!(
        ax,
        times / day,
        getproperty(default_timeseries, sym);
        label="Vopt = 10",
        linewidth=2,
    )
    lines!(
        ax, times / day, getproperty(vopt_timeseries, sym); label="Vopt = 5", linewidth=2
    )
    axislegend(ax; position=:rt)
end

#Save figure
save("predator_optimal_prey_size.png", fig; px_per_unit=1)

fig  # Display the figure
