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

# Define all plankton entities in a named tuple
plankton = (
    P1=(volume=1, opt_ratio=0, protection=0),
    P2=(volume=10, opt_ratio=0, protection=0),
    Z1=(volume=10, opt_ratio=10, protection=1),
    Z2=(volume=100, opt_ratio=10, protection=1),
)

# Define a function to create the palatability matrix
function create_palatability_matrix(plankton, volume_key, opt_ratio_key, protection_key)
    palatability_matrix = [
        dummy_emergent_palat(
            plankton[prey_name][volume_key],
            plankton[pred_name][volume_key],
            plankton[pred_name][opt_ratio_key],
            plankton[prey_name][protection_key],
        ) for pred_name in keys(plankton), prey_name in keys(plankton)
    ]
    return palatability_matrix
end

# Use the function to create the palatability matrix
palatability_matrix = create_palatability_matrix(plankton, :volume, :opt_ratio, :protection)

print(palatability_matrix) #should be [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 1 0.3 0.0 0.0; 0.3 1 0 0.0]
