function estimate_emergent_dictionary(plankton, growth_function, params)
    emergent_vector = Dict(
        name => growth_function(
            [plankton[name][key] for key in params]...,  # Use splatting to unpack the vector
        ) for name in keys(plankton)
    )
    return emergent_vector
end


using NamedArrays

# Define the dynamic palatability matrix creation function
function estimate_emergent_matrix(plankton, func, keys)
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
            arguments = [plankton[prey_name][k] for k in keys]  # dynamic prey keys
            arguments = vcat(arguments, [plankton[pred_name][k] for k in keys])  # dynamic predator keys

            # Calculate the palatability using the provided function and arguments
            palatability_values[i, j] = func(arguments...)
        end
    end

    # Create and return a NamedArray with the palatability values and names
    return NamedArray(palatability_values, names)
end

