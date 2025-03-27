# # [Zooplankton immunity] (@id immunity_example)

using Agate
using Agate.Constructors: NiPiZD
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units

const year = years = 365day
using Graphs, GraphPlot, Colors, ColorSchemes

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

# Function to map values to Viridis colors
function value_to_color(value, min_val, max_val)
    norm_value = (value - min_val) / (max_val - min_val + 1e-5)  # Avoid division by zero
    return get(ColorSchemes.viridis, norm_value)  # Use Viridis color scheme
end

function plot_graph_with_colored_edges(g, matrix, labels)
    edge_list = collect(edges(g))  # Extract edges as SimpleEdge{Int}
    edge_values = [matrix[e.dst, e.src] for e in edge_list]  # Extract edge weights from the TRANSPOSED matrix
    min_val, max_val = minimum(edge_values), maximum(edge_values)
    edge_colors = [value_to_color(v, min_val, max_val) for v in edge_values]    
    gplot(g,
        nodelabel=labels,  # Add node labels
        edgestrokec=edge_colors,  # Color edges based on weights
        edgelabel=edge_values, 
        layout=circular_layout
    )
end

# ## Plot default palatability matrices
g_palat = SimpleDiGraph(palat_matrix.array')
plot_graph_with_colored_edges(g_palat, palat_matrix.array, labels)

# ## Update palatability matrix

# Add new edges:
z1, z2 = findfirst(isequal("Z1"), labels), findfirst(isequal("Z2"), labels)
add_edge!(g_palat, z1, z2)
palat_matrix.array .= Graphs.adjacency_matrix(g_palat)' # ' to transpose
nothing #hide

# Plot:
plot_graph_with_colored_edges(g_palat, palat_matrix.array, labels)


#  ## Update assimilation matrix
# We also need to update the assimilation matrix. 
# We set the value for the new link to be the default 0.32.

g_assim = SimpleDiGraph(assim_matrix.array')
add_edge!(g_assim, z1, z2)
assim_matrix.array .= Graphs.adjacency_matrix(g_palat)' # ' to transpose
nothing #hide



# ## Instantiate models
new_model = NiPiZD.instantiate(N2P2ZD; 
    palatability_matrix=palat_matrix,
    assimilation_efficiency_matrix=assim_matrix)


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


# # Plotting
# using Plots

# fields = [:P1, :P2, :Z1, :Z2, :D, :N]
# titles = [
#     "Phytoplankton 1",
#     "Phytoplankton 2",
#     "Zooplankton 1",
#     "Zooplankton 2",
#     "Detritus",
#     "Nutrient",
# ]

# plots = [
#     plot(; title=titles[i], xlabel="Time (days)", ylabel=string(fields[i])) for
#     i in 1:length(fields)
# ]

# for (i, field) in enumerate(fields)
#     plot!(plots[i], timeseries_default[field]; label="immunity", color=:blue)
#     plot!(plots[i], timeseries_new[field]; label="no immunity", color=:green)
# end

# p = plot(plots...; layout=(2, 3), size=(900, 600))
# p
