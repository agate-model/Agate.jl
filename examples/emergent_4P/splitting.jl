# Define the function for dynamic parameter parsing
function parse_parameters(defined_parameters)
    # Initialize the resulting dictionary
    intermediate_parameters = Dict()

    # Iterate over each key in defined_parameters (e.g., "P", "Z", "cocco")
    for (key, params) in defined_parameters
        # Retrieve parameters specific to this key
        n = params["n"]
        min_volume = params["min_volume"]
        max_volume = params["max_volume"]

        # Get the splitting function from the dictionary and apply it to calculate volumes
        splitting_function = params["splitting"]
        volumes = splitting_function(min_volume, max_volume, n)

        # Create sub-dictionaries for each instance in this category
        for i in 1:n
            # Create a sub-dictionary for each "P1", "P2", etc.
            sub_dict = Dict()

            # Set the volume parameter for this entry
            sub_dict["volume"] = volumes[i]

            # Dynamically add all other parameters except the special ones (n, min_volume, max_volume, splitting)
            for (param_key, value) in params
                if !(param_key in ["n", "min_volume", "max_volume", "splitting"])
                    sub_dict[param_key] = value
                end
            end

            # Insert the sub-dictionary into intermediate_parameters with dynamic key names (e.g., "P1", "Z1")
            intermediate_parameters["$key$i"] = sub_dict
        end
    end

    return intermediate_parameters
end

# Define the splitting functions
function log_splitting(min_volume, max_volume, n)
    log_min = log10(min_volume)
    log_max = log10(max_volume)
    log_step = (log_max - log_min) / (n - 1)
    return [10^(log_min + i * log_step) for i in 0:(n - 1)]
end

function linear_splitting(min_volume, max_volume, n)
    linear_step = (max_volume - min_volume) / (n - 1)
    return [min_volume + i * linear_step for i in 0:(n - 1)]
end

# Example usage with embedded splitting functions
defined_parameters = Dict(
    "P" => Dict(
        "n" => 4,
        "min_volume" => 1,
        "max_volume" => 100,
        "splitting" => log_splitting,
        "protection" => 0,
        "growth_a" => 2.5,
    ),
    "Z" => Dict(
        "n" => 4,
        "min_volume" => 100,
        "max_volume" => 1000,
        "splitting" => linear_splitting,
        "protection" => 1,
        "growth_a" => 25,
    ),
)

# Generate intermediate parameters using embedded splitting functions
intermediate_parameters = parse_parameters(defined_parameters)

# Display the results
println("Intermediate Parameters:")
println(intermediate_parameters)
