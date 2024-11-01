using NamedArrays

# Define the dynamic palatability matrix creation function
function emergent_matrix(plankton, func, params)
    # Extract predator and prey names
    predator_names = collect(keys(plankton))
    prey_names = collect(keys(plankton))

    # Initialize a NamedArray to hold palatability values
    palatability_values = zeros(Float64, length(predator_names), length(prey_names))
    names = (predator=predator_names, prey=prey_names)

    # Populate the NamedArray with calculated values
    for (i, pred_name) in enumerate(predator_names)
        for (j, prey_name) in enumerate(prey_names)
            # Extract values dynamically for both predator and prey based on provided keys
            arguments = [plankton[prey_name][k] for k in params]  # dynamic prey keys
            arguments = vcat(arguments, [plankton[pred_name][k] for k in params])  # dynamic predator keys

            # Calculate the palatability using the provided function and arguments
            palatability_values[i, j] = func(arguments...)
        end
    end

    # Create and return a NamedArray with the palatability values and names
    return NamedArray(palatability_values, names)
end

# Function to generate a dictionary of emergent values for each species
function emergent_dict(plankton, func, params)
    emergent_values = Dict(
        name => func(
            [plankton[name][key] for key in params]...,  # Extract parameters dynamically
        ) for name in keys(plankton)
    )
    return emergent_values
end

# Main function that infers the type of output by running the model and checking output
function emergent_analysis(plankton, func, params)
    # Run both matrix and dictionary models, and check which output type is correct
    matrix_result = emergent_matrix(plankton, func, params)
    dict_result = emergent_dict(plankton, func, params)

    # Check which output is valid by type and dispatch the appropriate handling function
    if isa(matrix_result, NamedArray)
        handle_result(matrix_result)
    elseif isa(dict_result, Dict)
        handle_result(dict_result)
    else
        error("Unexpected output type: Expected either NamedArray or Dict.")
    end
end

# Handle NamedArray result (matrix logic)
function handle_result(result::NamedArray)
    println("Handling NamedArray (Matrix) result")
    return display(result)
end

# Handle Dict result (dictionary logic)
function handle_result(result::Dict)
    println("Handling Dict result")
    return display(result)
end

# Test with sample data
plankton = Dict(
    "species1" =>
        Dict("value" => 1.0, "pred" => 2.0, "opt_ratio" => 0.5, "protection" => 0.1),
    "species2" =>
        Dict("value" => 1.5, "pred" => 1.8, "opt_ratio" => 0.6, "protection" => 0.2),
)
keys = ["value", "pred", "opt_ratio", "protection"]
params = ["value", "pred", "opt_ratio", "protection"]

# Define sample functions for matrix and dictionary generation
palatability_func(x, y, opt_ratio, protection) = x * y * opt_ratio - protection
dict_func(args...) = sum(args)

# Call emergent_analysis for both cases, which should automatically infer correct handling
println("Testing with Matrix-generating function:")
emergent_analysis(plankton, palatability_func, keys)  # Should infer matrix result

println("\nTesting with Dict-generating function:")
emergent_analysis(plankton, dict_func, params)  # Should infer dictionary result
