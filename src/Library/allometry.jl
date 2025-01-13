module Allometry

export allometric_scaling_power,
    allometric_palatability_unimodal, allometric_palatability_unimodal_protection

"""
    allometric_scaling_power(a, b, d::Number)

Allometric scaling function using the power law for cell volume.

# Arguments
- `a`: scale
- `b`: exponent
- `d`: cell equivalent spherical diameter (ESD)
"""
function allometric_scaling_power(a::Number, b::Number, d::Number)
    V = (4 / 3) * π * (d / 2)^3
    return a * V^b
end

"""
    allometric_palatability_unimodal(prey_data, predator_data)

Calculates the unimodal allometric palatability of prey based on predator-prey diameters.

This function extracts `prey_diameter`, `predator_diameter`, `optimum_predator_prey_ratio`, 
and `specificity` from the provided dictionaries and calculates the palatability using the diameter-based formula.

Note that this formulation differs from the currently operational MITgcm-DARWIN model as it uses diameter instead of volumes and is structurally different. However, both formulations result in a unimodal response the width and optima are modulated by the optimum_predator_prey ratio and the specificity.

# Arguments
- `prey_data`: A dictionary containing prey-specific data:
  - `diameters`: Diameter of the prey.
- `predator_data`: A dictionary containing predator-specific data:
  - `diameters`: Diameter of the predator.
  - `optimum_predator_prey_ratio`: The optimal predator-prey diameter ratio for the predator.
  - `specificity`: A parameter controlling how sharply the palatability decreases away from the optimal ratio.

# Returns
- `palatability`: A number between 0 and 1 representing the palatability.
"""
function allometric_palatability_unimodal(prey_data::Dict, predator_data::Dict)
    prey_diameter = prey_data["diameters"]
    predator_diameter = predator_data["diameters"]
    predator_prey_optimum = predator_data["optimum_predator_prey_ratio"]
    predator_specificity = predator_data["specificity"]
    can_eat = predator_data["can_eat"]

    if can_eat == 0
        palatability = 0
    elseif can_eat == 1
        predator_prey_ratio = prey_diameter / predator_diameter
        palatability =
            1 / (1 + (predator_prey_ratio - predator_prey_optimum)^2)^predator_specificity
    end
    return palatability
end

"""
    allometric_palatability_unimodal_protection(prey_data, predator_data)

Calculates the unimodal allometric palatability of prey, accounting for additional prey protection mechanisms.

The function uses a modified unimodal relationship defined by:
`palatability = prey_protection / (1 + (predator_prey_ratio - predator_prey_optimum)^2)^predator_specificity`

# Arguments
- `prey_data`: A dictionary containing prey-specific data:
  - `diameters`: Diameter of the prey.
  - `protection`: A scaling factor between 0 and 1 representing additional protection mechanisms of the prey.
- `predator_data`: A dictionary containing predator-specific data:
  - `diameters`: Diameter of the predator.
  - `optimum_predator_prey_ratio`: The optimal predator-prey diameter ratio for the predator.
  - `specificity`: A parameter controlling how sharply the palatability decreases away from the optimal ratio.

# Returns
- `palatability`: A number between 0 and `prey_protection` representing the palatability.
"""
function allometric_palatability_unimodal_protection(prey_data::Dict, predator_data::Dict)
    prey_diameter = prey_data["diameters"]
    predator_diameter = predator_data["diameters"]
    predator_prey_optimum = predator_data["optimum_predator_prey_ratio"]
    predator_specificity = predator_data["specificity"]
    prey_protection = prey_data["protection"]
    can_eat = predator_data["can_eat"]

    if can_eat == 0
        palatability = 0
    elseif can_eat == 1
        predator_prey_ratio = predator_diameter / prey_diameter
        palatability =
            (1 - prey_protection) /
            (1 + (predator_prey_ratio - predator_prey_optimum)^2)^predator_specificity
    end

    return palatability
end

end
