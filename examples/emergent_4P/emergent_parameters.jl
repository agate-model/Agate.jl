using NamedArrays

# Define a dummy function for calculating growth
function dummy_emergent_growth(growth_a, growth_b, volume)
    rate = 0.0
    if growth_a == 0
        return rate
    else
        if volume == 1
            rate = 7.190e-6
        elseif volume == 10
            rate = 2.216e-5
        end
    end

    return rate
end

# Define a dummy function for calculating palatability
function dummy_emergent_palat(
    prey_volume, predator_volume, optimum_predator_prey_ratio, protection
)
    ratio = predator_volume / prey_volume

    if optimum_predator_prey_ratio == 0
        palat = 0.0
    elseif ratio == optimum_predator_prey_ratio
        palat = 1 * (1 - protection)
    else
        palat = 0.3 * (1 - protection)
    end
    return palat
end

function estimate_emergent_vectors(plankton, growth_function, params)
    emergent_vector = Dict(
        name => growth_function(
            [plankton[name][key] for key in params]...  # Use splatting to unpack the vector
        ) for name in keys(plankton)
    )
    return emergent_vector
end

# Define a function to create the palatability matrix using NamedArrays
function create_palatability_matrix(plankton, volume_key, opt_ratio_key, protection_key)
    # Extract predator and prey names
    predator_names = collect(keys(plankton))
    prey_names = collect(keys(plankton))

    # Initialize a NamedArray to hold palatability values
    palatability_values = zeros(Float64, length(predator_names), length(prey_names))
    names = (predator=predator_names, prey=prey_names)

    # Populate the NamedArray with calculated values
    for (i, pred_name) in enumerate(predator_names)
        for (j, prey_name) in enumerate(prey_names)
            # Calculate the palatability using the provided function
            palatability_values[i, j] = dummy_emergent_palat(
                plankton[prey_name][volume_key],
                plankton[pred_name][volume_key],
                plankton[pred_name][opt_ratio_key],
                plankton[prey_name][protection_key],
            )
        end
    end

    # Create a NamedArray with the palatability values and names
    palatability_matrix = NamedArray(palatability_values, names)
    
    return palatability_matrix
end

# Example plankton data
plankton = Dict(
    "P1" => Dict(
        "volume" => 1,
        "growth_a" => 1,
        "growth_b" => 1,
        "protection" => 0,
        "opt_ratio" => 0,
    ),
    "P2" => Dict(
        "volume" => 10,
        "growth_a" => 1,
        "growth_b" => 1,
        "protection" => 0,
        "opt_ratio" => 0,
    ),
    "Z1" => Dict(
        "volume" => 10,
        "growth_a" => 0,
        "growth_b" => 0,
        "protection" => 1,
        "opt_ratio" => 10,
    ),
    "Z2" => Dict(
        "volume" => 100,
        "growth_a" => 0,
        "growth_b" => 0,
        "protection" => 1,
        "opt_ratio" => 10,
    ),
)

# Using the updated function with completely dynamic arguments
growth_rates = estimate_emergent_vectors(
    plankton,
    dummy_emergent_growth,  # Function to call
    ["growth_a", "growth_b", "volume"],  # Use strings for parameter keys
)

println("growth_rates: ", growth_rates)

# Create the palatability matrix
palatability_matrix = create_palatability_matrix(
    plankton, "volume", "opt_ratio", "protection"
)
