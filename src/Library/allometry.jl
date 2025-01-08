module Allometry

"""
    allometric_scaling_power(a, b, V::Number)

Allometric scaling function using the power law for cell volume.

# Arguments
- `a`: scale
- `b`: exponent
- `V`: cell volume
"""
function allometric_scaling_power(a::Number, b::Number, V::Number)
    return a * V^b
end

"""
    allometric_scaling_power(a, b, d::Number)

Allometric scaling function using the power law for cell diameter. Converts
diameter to volume assuming a spherical shape.

Note that if diameter is passed instead of volume this should be done explicitly:
`x = (2.19, -0.16, d=10)`

# Arguments
- `a`: scale
- `b`: exponent
- `d`: cell diameter
"""
function allometric_scaling_power(a::Number, b::Number, d::Number)
    V = (4 / 3) * π * (d / 2)^3
    return allometric_scaling_power(a, b, V)
end

"""
    allometric_palatability_unimodal(prey_volume, predator_volume, predator_prey_optimum, predator_specificity)

Calculates the unimodal allometric palatability of prey based on the predator-prey volume ratio.

The function uses a unimodal relationship defined by:
`palatability = 1 / (1 + (predator_prey_ratio - predator_prey_optimum)^2)^predator_specificity`

# Arguments
- `prey_volume`: Volume of the prey (assumed to be spherical if diameter is provided).
- `predator_volume`: Volume of the predator (assumed to be spherical if diameter is provided).
- `predator_prey_optimum`: The optimal predator-prey ratio for the predator.
- `predator_specificity`: A parameter controlling how sharply the palatability decreases away from the optimal ratio.

# Returns
- `palatability`: A number between 0 and 1 representing the palatability.
"""
function allometric_palatability_unimodal(
    prey_volume::Number,
    predator_volume::Number,
    predator_prey_optimum::Number,
    predator_specificity::Number,
)
    predator_prey_ratio = prey_volume / predator_volume
    palatability =
        1 / (1 + (predator_prey_ratio - predator_prey_optimum)^2)^predator_specificity
    return palatability
end

"""
    allometric_palatability_unimodal(prey_diameter, predator_diameter, predator_prey_optimum, predator_specificity)

Calculates the unimodal allometric palatability of prey based on the predator-prey diameter ratio.

This method converts prey and predator diameters into their respective volumes (assuming spherical geometry) and uses the volume-based formula.

# Arguments
- `prey_diameter`: Diameter of the prey (assumed to be spherical).
- `predator_diameter`: Diameter of the predator (assumed to be spherical).
- `predator_prey_optimum`: The optimal predator-prey ratio for the predator.
- `predator_specificity`: A parameter controlling how sharply the palatability decreases away from the optimal ratio.

# Returns
- `palatability`: A number between 0 and 1 representing the palatability.
"""
function allometric_palatability_unimodal(
    prey_diameter::Number,
    predator_diameter::Number,
    predator_prey_optimum::Number,
    predator_specificity::Number,
)
    prey_volume = (4 / 3) * π * (prey_diameter / 2)^3
    predator_volume = (4 / 3) * π * (predator_diameter / 2)^3
    return allometric_palatability_unimodal(
        prey_volume, predator_volume, predator_prey_optimum, predator_specificity
    )
end

"""
    allometric_palatability_unimodal_protection(prey_volume, predator_volume, predator_prey_optimum, predator_specificity, prey_protection)

Calculates the unimodal allometric palatability of prey, accounting for additional prey protection mechanisms.

The function uses a modified unimodal relationship defined by:
`palatability = prey_protection / (1 + (predator_prey_ratio - predator_prey_optimum)^2)^predator_specificity`

# Arguments
- `prey_volume`: Volume of the prey (assumed to be spherical if diameter is provided).
- `predator_volume`: Volume of the predator (assumed to be spherical if diameter is provided).
- `predator_prey_optimum`: The optimal predator-prey ratio for the predator.
- `predator_specificity`: A parameter controlling how sharply the palatability decreases away from the optimal ratio.
- `prey_protection`: A scaling factor between 0 and 1 representing additional protection mechanisms of the prey.

# Returns
- `palatability`: A number between 0 and `prey_protection` representing the palatability.
"""
function allometric_palatability_unimodal_protection(
    prey_volume::Number,
    predator_volume::Number,
    predator_prey_optimum::Number,
    predator_specificity::Number,
    prey_protection::Number,
)
    return prey_protection / allometric_palatability_unimodal(
        prey_volume, predator_volume, predator_prey_optimum, predator_specificity
    )
end

"""
    allometric_palatability_unimodal_protection(prey_diameter, predator_diameter, predator_prey_optimum, predator_specificity, prey_protection)

Calculates the unimodal allometric palatability of prey, accounting for additional prey protection mechanisms, based on diameters.

This method converts prey and predator diameters into their respective volumes (assuming spherical geometry) and uses the volume-based formula.

# Arguments
- `prey_diameter`: Diameter of the prey (assumed to be spherical).
- `predator_diameter`: Diameter of the predator (assumed to be spherical).
- `predator_prey_optimum`: The optimal predator-prey ratio for the predator.
- `predator_specificity`: A parameter controlling how sharply the palatability decreases away from the optimal ratio.
- `prey_protection`: A scaling factor between 0 and 1 representing additional protection mechanisms of the prey.

# Returns
- `palatability`: A number between 0 and `prey_protection` representing the palatability.
"""
function allometric_palatability_unimodal_protection(
    prey_diameter::Number,
    predator_diameter::Number,
    predator_prey_optimum::Number,
    predator_specificity::Number,
    prey_protection::Number,
)
    prey_volume = (4 / 3) * π * (prey_diameter / 2)^3
    predator_volume = (4 / 3) * π * (predator_diameter / 2)^3
    return allometric_palatability_unimodal_protection(
        prey_volume,
        predator_volume,
        predator_prey_optimum,
        predator_specificity,
        prey_protection,
    )
end

end
