# # [Zooplankton immunity] (@id immunity_example)

using Agate
using Agate.Constructors: NiPiZD
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units
using Graphs, GraphPlot

const year = years = 365day


# Construct the ecosystem model
N2P2ZD = NiPiZD.construct()
default_model = N2P2ZD()
nothing #hide

# Extract palatability matrix
palat_matrix = default_model.palatability_matrix
assim_matrix = default_model.assimilation_efficiency_matrix

println(assim_matrix)
labels = collect(keys(default_model.assimilation_efficiency_matrix.dicts[1]))
nothing #hide

# Turn matrix into graph and plot
# !!! info
#     In Agate.jl we follow the ecological predator-prey matrix convention where predators = rows, and prey = columns.
#     However, Graphs.jl follows the graph theory convention where the sources (prey) = rows, and the targets (predators) = columns.
#     To plot the graphs, and to add nodes, we thus need to transpose (`'`) the matrices before plotting and before instantiating the new model.


g = SimpleDiGraph(palat_matrix.array') # ' to transpose
gplot(g, nodelabel=labels, layout=circular_layout)

# Add predation of Z2 on Z1
z1, z2 = findfirst(isequal("Z1"), labels), findfirst(isequal("Z2"), labels)
add_edge!(g, z1, z2)
nothing #hide

# Plot the new graph
gplot(g, nodelabel=labels, layout=circular_layout)

# Update palat_matrix with the new adjacency matrix
palat_matrix.array .= Graphs.adjacency_matrix(g)' # ' to transpose
nothing #hide

# Instantiate the updated model
new_model = NiPiZD.instantiate(N2P2ZD; palatability_matrix=palat_matrix)


bgc_model_default = Biogeochemistry(
    default_model;
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()),
)
full_model_default = BoxModel(;
    biogeochemistry=bgc_model_default
)
set!(full_model_default; P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)



bgc_model_new = Biogeochemistry(
    new_model;
    light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid()),
)
full_model_new = BoxModel(;
    biogeochemistry=bgc_model_new
)
set!(full_model_new; P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, N=7.0, D=1)


# Simulate models

function run_simulation!(model, filename)
    simulation = Simulation(model; Δt=10minutes, stop_time=1years)
    simulation.output_writers[:fields] = JLD2OutputWriter(
        model,
        model.fields;
        filename=filename,
        schedule=TimeInterval(1day),
        overwrite_existing=true,
    )
    run!(simulation)
    return NamedTuple{keys(model.fields)}(
        FieldTimeSeries(filename, "$field")[1, 1, 1, :] for field in keys(model.fields)
    )
end

filename_default = "box_n2p2zd_default.jld2"
filename_new = "box_n2p2zd_new.jld2"

timeseries_default = run_simulation!(
    full_model_default, filename_default
)
timeseries_new = run_simulation!(
    full_model_new, filename_new
)


# ==================================================
# Plotting
# ==================================================
using Plots

fields = [:P1, :P2, :Z1, :Z2, :D, :N]
titles = [
    "Phytoplankton 1",
    "Phytoplankton 2",
    "Zooplankton 1",
    "Zooplankton 2",
    "Detritus",
    "Nutrient",
]

plots = [
    plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
    i in 1:length(fields)
]

for (i, field) in enumerate(fields)
    plot!(plots[i], timeseries_default[field]; label="immunity", color=:blue)
    plot!(plots[i], timeseries_new[field]; label="no immunity", color=:green)
end

p = plot(plots...; layout=(2, 3), size=(900, 600))
p
