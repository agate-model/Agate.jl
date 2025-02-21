"""
Functions related to zooplankton predation.
"""

module Predation

export holling_type_2,
    predation_loss_idealized,
    predation_gain_idealized,
    predation_assimilation_loss_idealized,
    predation_loss_preferential,
    predation_gain_preferential,
    predation_assimilation_loss_preferential,
    assimilation_efficiency_emergent_binary

"""
Holling's "type II" functional response as describe in Holling 1959.
The function is similar to the Monod equation and Michaelis-Menten equation of for enzyme kinetics.
The formulation is characterized by decelerating predation as prey concentrations increase.

# Arguments
- `R`: prey density
- `k`: prey density at which predation is half it's maximum rate
"""
function holling_type_2(R::Real, k::Real)
    return R / (k + R)
end

#idealized predation (OceanBioME NPZD model)

"""
Estimates the loss rate of P (prey), to Z (predator).
In this formulation predator-prey interactions are modulated both by their density (Holling type 2)
and the prey-predator palatability.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: prey density at which predation is half it's maximum rate
"""
function predation_loss_idealized(P, Z, gₘₐₓ, kₚ)
    return gₘₐₓ * holling_type_2(P^2, kₚ^2) * Z
end

"""
Estimates the gain rate of Z (predator) feeding on P (prey).
In this formulation predation rate is multiplied by an assimilation efficiency (β) which
represents 'sloppy feeding'.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: grazing/holling half saturation
"""
function predation_gain_idealized(P, Z, β, gₘₐₓ, kₚ)
    return β * predation_loss_idealized(P, Z, gₘₐₓ, kₚ)
end

"""
Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

Note that this differs from the predation_gain_idealized as this function represents the transfer of biomass from the prey to the environment
rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing/holling half saturation
"""
function predation_assimilation_loss_idealized(P, Z, β, gₘₐₓ, kₚ)
    return (1 - β) * predation_loss_idealized(P, Z, gₘₐₓ, kₚ)
end

#preferential predation (intermediate complexity model)

"""
Estimates the loss rate of P (prey), to Z (predator).
In this formulation predator-prey interactions are modulated both by their density (Holling type 2)
and the prey-predator palatability.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: prey density at which predation is half it's maximum rate
- `palatability`: the likelihood at which the predator feeds on the prey
"""
function predation_loss_preferential(P, Z, gₘₐₓ, kₚ, palatability)
    return gₘₐₓ * palatability * holling_type_2(P, kₚ) * Z
end

"""
Estimates the gain rate of Z (predator) feeding on P (prey).
In this formulation predation rate is multiplied by an assimilation efficiency (β) which
represents 'sloppy feeding'.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: grazing/holling half saturation
- `palatability`: the likelihood at which the predator feeds on the prey
"""
function predation_gain_preferential(P, Z, β, gₘₐₓ, kₚ, palatability)
    return β * predation_loss_preferential(P, Z, gₘₐₓ, kₚ, palatability)
end

"""
Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

Note that this differs from the predation_gain_preferential as this function represents the transfer of biomass from the prey to the environment
rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing/holling half saturation
- `palatability`: the likelihood at which the predator feeds on the prey
"""
function predation_assimilation_loss_preferential(P, Z, β, gₘₐₓ, kₚ, palatability)
    return (1 - β) * predation_loss_preferential(P, Z, gₘₐₓ, kₚ, palatability)
end

"""
    assimilation_efficiency_emergent_binary(prey_data, predator_data)

Determines the assimilation efficiency of a predator consuming prey, based on binary conditions of edibility.

The function evaluates whether the predator can eat the prey and whether the prey can be consumed, and assigns the assimilation efficiency accordingly.

# Arguments
- `prey_data`: A dictionary containing prey-specific data:
  - `can_be_eaten`: A binary value (1 or 0) indicating if the prey can be consumed by the predator.
- `predator_data`: A dictionary containing predator-specific data:
  - `can_eat`: A binary value (1 or 0) indicating if the predator can consume prey.
  - `assimilation_efficiency`: The efficiency with which the predator assimilates nutrients from the prey if the conditions are met.

# Returns
- `assimilation_efficiency`:
  - If `can_eat` is 1 and `can_be_eaten` is 1, returns the predator's `assimilation_efficiency`.
  - Otherwise, returns 0.
"""
function assimilation_efficiency_emergent_binary(prey_data, predator_data)
    if predator_data["can_eat"] == 1 && prey_data["can_be_eaten"] == 1
        assimilation_efficiency = predator_data["assimilation_efficiency"]
    elseif predator_data["can_eat"] == 1 && prey_data["can_be_eaten"] == 0
        assimilation_efficiency = 0
    elseif predator_data["can_eat"] == 0
        assimilation_efficiency = 0
    end
    return assimilation_efficiency
end

end # module
