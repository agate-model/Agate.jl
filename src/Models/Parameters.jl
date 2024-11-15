module Parameters

using NamedArrays

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
Generate volumes....

Key (e.g., "P", "Z", "cocco")
"""
function split_size_parameters(
    plankton_name, n_plankton, min_volume, max_volume, splitting_function=linear_splitting
)
    volumes = splitting_function(min_volume, max_volume, n_plankton)

    # Initialize the resulting dictionary
    intermediate_parameters = Dict()
    # Create sub-dictionaries for each instance in this category
    for i in 1:n_plankton
        # Create a sub-dictionary for each "P1", "P2", etc.
        sub_dict = Dict()

        # Set the volume parameter for this entry
        sub_dict["volume"] = volumes[i]

        intermediate_parameters["$plantkon_name$i"] = sub_dict
    end

    return intermediate_parameters
end

"""
To use with:
    - maximum_growth_rate
    - maximum_predation_rate
    - nitrogen_half_saturation

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

function palatability(plankton::Dict, func::Function, key_list::Vector{String})
    # Extract predator and prey names
    predator_names = collect(keys(plankton))
    prey_names = collect(keys(plankton))

    # Initialize a NamedArray to hold palatability values
    palatability_values = zeros(Float64, length(predator_names), length(prey_names))
    names = (predator=predator_names, prey=prey_names)

    # Populate the NamedArray with calculated values
    for (i, pred_name) in enumerate(predator_names)
        for (j, prey_name) in enumerate(prey_names)
            prey_data = plankton[prey_name]
            predator_data = plankton[pred_name]

            # Pass prey and predator data dictionaries to the function with dynamic keys
            palatability_values[i, j] = func(
                prey_data,
                predator_data;
                prey_volume_key=key_list[1],
                predator_volume_key=key_list[1],
                optimum_ratio_key=key_list[2],
                protection_key=key_list[3],
            )
        end
    end

    # Create and return a NamedArray with the palatability values and names
    return NamedArray(palatability_values, names)
end

end #module