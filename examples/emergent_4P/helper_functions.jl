using NamedArrays

"""

#Example usage:

plankton = Dict(
    "species1" =>
        Dict("volume" => 1.0,"volume_a" => 1.0,"volume_b" => 1.0, "pred" => 2.0),
    "species2" =>
        Dict("volume" => 1.5,"volume_a" => 1.0,"volume_b" => 1.0, "pred" => 1.8),
)

result = emergent_1D_array(plankton, dummy_emergent_predation_rate, ["volume_a", "volume_b", "volume"])
"""
function emergent_1D_array(plankton::Dict, func::Function, params::Vector{String})
    # Get species names
    species_names = collect(keys(plankton))

    # Calculate emergent values for each species
    emergent_values = [
        func([plankton[name][key] for key in params]...) for name in species_names
    ]

    # Create a NamedArray with species names as row labels
    emergent_array = NamedArray(emergent_values, (species_names,), ("Species",))

    return emergent_array
end

"""
#Example usage:

# Sample plankton data with required keys
plankton = Dict(
    "species1" =>
        Dict("volume" => 1.0, "pred" => 2.0, "optimum_predator_prey_ratio" => 10, "protection" => 0),
    "species2" =>
        Dict("volume" => 1.5, "pred" => 1.8, "optimum_predator_prey_ratio" => 10, "protection" => 0),
)

# Specify keys for prey volume, predator volume, optimum ratio, and protection
result = emergent_matrix(plankton, dummy_emergent_palat, ["volume", "optimum_predator_prey_ratio", "protection"])

"""
function emergent_2D_array(plankton::Dict, func::Function, key_list::Vector{String})
    # Extract predator and prey names
    predator_names = collect(keys(plankton))
    prey_names = collect(keys(plankton))

    # Initialize a NamedArray to hold values
    values = zeros(Float64, length(predator_names), length(prey_names))
    names = (predator=predator_names, prey=prey_names)

    # Populate the NamedArray with calculated values
    for (i, pred_name) in enumerate(predator_names)
        for (j, prey_name) in enumerate(prey_names)
            prey_data = plankton[prey_name]
            predator_data = plankton[pred_name]
            # println("prey_data", predator_data)
            # Pass prey and predator data dictionaries to the function
            values[i, j] = func(prey_data, predator_data;)
        end
    end

    # Create and return a NamedArray with the values and names
    return NamedArray(values, names)
end

"""
#Example:

emergent_analysis(
    plankton, dummy_emergent_palat, ["volume", "optimum_predator_prey_ratio", "protection"]
)

"""
function emergent_analysis(plankton::Dict, func::Function, params::Vector{String})
    result = nothing
    # println(func)
    try
        result = emergent_1D_array(plankton, func, params)
        if ndims(result_1D) == 1
            println("1D")
        else
            error("Unexpected output type")
        end
    catch
    end
    try
        result = emergent_2D_array(plankton, func, params)

        if ndims(result_2D) == 2
            println("2D")
        else
            error("Unexpected output type")
        end
    catch
    end

    return result
end

"""
A function which takes a dictionary containing n, min_size, max_size, a splitting function,
and other biogeochemical parameters.

#Example:

defined_parameters = Dict(
    "P" => Dict(
        "n" => 4,
        "min_volume" => 1,
        "max_volume" => 100,
        "splitting" => log_splitting, #function should be defined
        "protection" => 0,  #aux variables
        "growth_a" => 2.5,  #aux variables
    ),
    "Z" => Dict(
        "n" => 4,
        "min_volume" => 100,
        "max_volume" => 1000,
        "splitting" => linear_splitting, #function should be defined
        "protection" => 1, #aux variables
        "growth_a" => 25,  #aux variables
    ),
)

intermediate_parameters = split_size_parameters(defined_parameters)

"""
function split_size_parameters(defined_parameters::Dict)
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

"""
Log splitting function to generate a set of volumes based on min_volume, max_volume and n
"""
function log_splitting(min_volume::Real, max_volume::Real, n::Int)
    log_min = log10(min_volume)
    log_max = log10(max_volume)
    log_step = (log_max - log_min) / (n - 1)
    return [10^(log_min + i * log_step) for i in 0:(n - 1)]
end

"""
Linear splitting function to generate a set of volumes based on min_volume, max_volume and n
"""
function linear_splitting(min_volume::Real, max_volume::Real, n::Int)
    linear_step = (max_volume - min_volume) / (n - 1)
    return [min_volume + i * linear_step for i in 0:(n - 1)]
end

"""
Function to do final step of emergent analysis
"""
function analyze_and_merge(emergent_functions::Dict, intermediate_parameters::Dict)
    # Dictionary to store results
    emergent_parameters = Dict()

    # First, handle parameters defined in `emergent_functions`
    for (func_name, (func, param_names)) in emergent_functions
        println("Analyzing function: $func_name")

        # Run the analysis with the specified function and store the result
        result = emergent_analysis(intermediate_parameters, func, param_names)
        println(result)
        emergent_parameters[func_name] = result
    end

    # Next, copy parameters directly for those not in `emergent_functions`
    for param_name in keys(intermediate_parameters["P1"])  # Assumes all sub-dicts have the same keys
        if !haskey(emergent_functions, param_name)
            println("Copying parameter directly: $param_name")

            # Collect values across species in an array
            values = [
                intermediate_parameters[species][param_name] for
                species in keys(intermediate_parameters)
            ]
            species_names = collect(keys(intermediate_parameters))

            # Store as a NamedArray
            emergent_parameters[param_name] = NamedArray(
                values, (species_names,), ("Species",)
            )
        end
    end

    # Return the full emergent_parameters dictionary
    return emergent_parameters
end
