using Graphs
using GraphPlot
using Colors
using ColorSchemes

# Function to map values to Viridis colors
function value_to_color(value, min_val, max_val)
    norm_value = (value - min_val) / (max_val - min_val + 1e-5)  # Avoid division by zero
    return get(ColorSchemes.viridis, norm_value)  # Use Viridis color scheme
end

function plot_graph_with_colored_edges(g, matrix, labels)
    edge_list = collect(edges(g))  # Extract edges as SimpleEdge{Int}
    edge_values = [matrix[e.src, e.dst] for e in edge_list]  # Extract edge weights

    min_val, max_val = minimum(edge_values), maximum(edge_values)
    edge_colors = [value_to_color(v, min_val, max_val) for v in edge_values]

    # Increase node size and scale text size
    node_size = 0.2  # Larger node size (relative to default)
    text_size = 0.8 * node_size  # Text size is 0.8 * node size
    edgelabel = edge_values

    # Plot the graph using GraphPlot
    gplot(g,
        nodelabel=labels,  # Add node labels
        nodelabelsize=text_size,  # Scale text size
        edgestrokec=edge_colors,  # Color edges based on weights
        edgelabel=edgelabel, 
        layout=circular_layout
    )
end

# Example graph and matrix
g = SimpleDiGraph(3)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 3, 1)

matrix = [0 0.5 0; 0 0 0.8; 0.3 0 0]  # Example weight matrix
labels = ["A", "B", "C"]

# Plot graph with colored edges
plot_graph_with_colored_edges(g, matrix, labels)