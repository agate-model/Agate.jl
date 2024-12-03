"""
Functions related to zooplankton predation.
"""

module Predation

export holling_type_2,
    idealized_predation_loss,
    idealized_predation_gain,
    idealized_predation_assimilation_loss

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
function idealized_predation_loss(P, Z, gₘₐₓ, kₚ)
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
function idealized_predation_gain(P, Z, β, gₘₐₓ, kₚ)
    return β * idealized_predation_loss(P, Z, gₘₐₓ, kₚ)
end

"""
Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

Note that this differs from the idealized_predation_gain as this function represents the transfer of biomass from the prey to the environment
rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing/holling half saturation
"""
function idealized_predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ)
    return (1 - β) * idealized_predation_loss(P, Z, gₘₐₓ, kₚ)
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
function preferential_predation_loss(P, Z, gₘₐₓ, kₚ, palatability)
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
function preferential_predation_gain(P, Z, β, gₘₐₓ, kₚ, palatability)
    return β * preferential_predation_loss(P, Z, gₘₐₓ, kₚ, palatability)
end

"""
Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

Note that this differs from the preferential_predation_gain as this function represents the transfer of biomass from the prey to the environment
rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing/holling half saturation
- `palatability`: ...
"""
function preferential_predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ, palatability)
    return (1 - β) * preferential_predation_loss(P, Z, gₘₐₓ, kₚ, palatability)
end

end # module
