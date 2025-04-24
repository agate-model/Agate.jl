module Allometry

export allometric_scaling_power,
    allometric_palatability_unimodal, allometric_palatability_unimodal_protection

"""
    allometric_scaling_power(a::Number, b::Number, diameter::Number)

Allometric scaling function using the power law for cell volume.

!!! formulation
    ``a````V````ᵇ``

    where:
    - ``V`` = (4 / 3) * π * (``d`` / 2)³
    - ``a`` = scale
    - ``b`` = exponent
    - ``d`` = cell equivalent spherical diameter (ESD)

# Arguments
- `a`: scale
- `b`: exponent
- `diameter`: cell equivalent spherical diameter (ESD)
"""
function allometric_scaling_power(a::Number, b::Number, diameter::Number)
    V = (4 / 3) * π * (diameter / 2)^3
    return a * V^b
end

"""
    allometric_palatability_unimodal(prey_data::Dict, predator_data::Dict)

Calculates the unimodal allometric palatability of prey based on predator-prey diameters.

!!! formulation
    0 if ``f`` = 0

    1 / (1 + (``d_{ratio}``- ``d_{opt}``)²)^σ   otherwise
    
    where:
    - ``f`` = binary ability of predator to eat prey
    - ``d_{ratio}`` = ratio between predator and prey diameters
    - ``d_{opt}`` = optimum ratio between predator and prey diameter
    - σ = how sharply the palatability decreases away from the optimal ratio.

!!! info   
    This formulation differs from the operational MITgcm-DARWIN model as it is is structurally different and diameters are used instead of volumes.
    However, both formulations result in a unimodal response where the width and optima are modulated by the optimum-predator-prey ratio and the specificity.

# Arguments
- `prey_data`: A dictionary containing prey-specific data:
  - `diameters`: Diameter of the prey.
- `predator_data`: A dictionary containing predator-specific data:
  - `can_eat`: A binary value (1 or 0) indicating if the predator can consume prey. If this is set to 0, palatability is set to 0.
  - `diameters`: Diameter of the predator.
  - `optimum_predator_prey_ratio`: The optimal predator-prey diameter ratio for the predator.
  - `specificity`: A parameter controlling how sharply the palatability decreases away from the optimal ratio.

# Returns
- `palatability`: A number between 0 and 1 representing the palatability.
"""
function allometric_palatability_unimodal(prey_data::Dict, predator_data::Dict)
    can_eat = predator_data["can_eat"]

    if can_eat == 0
        palatability = 0
    elseif can_eat == 1
        prey_diameter = prey_data["diameters"]
        predator_diameter = predator_data["diameters"]
        predator_prey_optimum = predator_data["optimum_predator_prey_ratio"]
        predator_specificity = predator_data["specificity"]
        predator_prey_ratio = prey_diameter / predator_diameter
        palatability =
            1 / (1 + (predator_prey_ratio - predator_prey_optimum)^2)^predator_specificity
    end
    return palatability
end

"""
    allometric_palatability_unimodal_protection(prey_data::Dict, predator_data::Dict)

Calculates the unimodal allometric palatability of prey, accounting for additional prey protection mechanisms.


!!! formulation
    0 if ``f`` = 0

    (1 - η) / (1 + (``d_{ratio}``- ``d_{opt}``)^2)^σ   otherwise
    
    where:
    - ``f`` = binary ability of predator to eat prey
    - η = prey-protection
    - ``d_{ratio}`` = ratio between predator and prey diameters
    - ``d_{opt}`` = optimum ratio between predator and prey diameter
    - σ = how sharply the palatability decreases away from the optimal ratio.


# Arguments
- `prey_data`: A dictionary containing prey-specific data:
  - `diameters`: Diameter of the prey.
  - `protection`: A scaling factor between 0 and 1 representing additional protection mechanisms of the prey.
- `predator_data`: A dictionary containing predator-specific data:
  - `can_eat`: A binary value (1 or 0) indicating if the predator can consume prey. If this is set to 0, palatability is set to 0.
  - `diameters`: Diameter of the predator.
  - `optimum_predator_prey_ratio`: The optimal predator-prey diameter ratio for the predator.
  - `specificity`: A parameter controlling how sharply the palatability decreases away from the optimal ratio.

# Returns
- `palatability`: A number between 0 and `prey_protection` representing the palatability.
"""
function allometric_palatability_unimodal_protection(prey_data::Dict, predator_data::Dict)
    can_eat = predator_data["can_eat"]

    if can_eat == 0
        palatability = 0
    elseif can_eat == 1
        prey_diameter = prey_data["diameters"]
        predator_diameter = predator_data["diameters"]
        predator_prey_optimum = predator_data["optimum_predator_prey_ratio"]
        predator_specificity = predator_data["specificity"]
        prey_protection = prey_data["protection"]

        predator_prey_ratio = predator_diameter / prey_diameter
        palatability =
            (1 - prey_protection) /
            (1 + (predator_prey_ratio - predator_prey_optimum)^2)^predator_specificity
    end

    return palatability
end

end
