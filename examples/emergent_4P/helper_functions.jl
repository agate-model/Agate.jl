using NamedArrays

```

#Example usage:

plankton = Dict(
    "species1" =>
        Dict("volume" => 1.0,"volume_a" => 1.0,"volume_b" => 1.0, "pred" => 2.0),
    "species2" =>
        Dict("volume" => 1.5,"volume_a" => 1.0,"volume_b" => 1.0, "pred" => 1.8),
)

result = emergent_1D_array(plankton, dummy_emergent_predation_rate, ["volume_a", "volume_b", "volume"])

```
function emergent_1D_array(plankton, func, params)
    # Get species names
    species_names = collect(keys(plankton))
    
    # Calculate emergent values for each species
    emergent_values = [func([plankton[name][key] for key in params]...) for name in species_names]
    
    # Create a NamedArray with species names as row labels
    emergent_array = NamedArray(
        emergent_values,
        (species_names,),
        ("Species",)
    )
    
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
function emergent_2D_array(plankton, func, key_list)
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
            palatability_values[i, j] = func(prey_data, predator_data;
                                             prey_volume_key=key_list[1],
                                             predator_volume_key=key_list[1],
                                             optimum_ratio_key=key_list[2],
                                             protection_key=key_list[3])
        end
    end

    # Create and return a NamedArray with the palatability values and names
    return NamedArray(palatability_values, names)
end
