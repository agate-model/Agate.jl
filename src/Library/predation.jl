"""
Functions related to zooplankton predation.
"""

module Predation

export holling_type_2,
    idealized_predation_loss,
    idealized_predation_gain,
    idealized_predation_assimilation_loss

"""
??

# Arguments
- `R`:
- `k`:
"""
function holling_type_2(R::Real, k::Real)
    return R / (k + R)
end

"""
Estimates the loss rate of P (prey), to Z (predator).

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: grazing/holling half saturation
- `palatability`: ...
"""
function idealized_predation_loss(P, Z, gₘₐₓ, kₚ, palatability)
    return gₘₐₓ * palatability * holling_type_2(P, kₚ) * Z
end

"""
Estimates the gain rate of Z (predator) feeding on P (prey).

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: grazing/holling half saturation
- `palatability`: ...
"""
function idealized_predation_gain(P, Z, β, gₘₐₓ, kₚ, palatability)
    return β * idealized_predation_loss(P, Z, gₘₐₓ, kₚ, palatability)
end

"""
Estimates the rate at which plankton predation gain is lost due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing/holling half saturation
- `palatability`: ...
"""
function idealized_predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ, palatability)
    return (1 - β) * idealized_predation_loss(P, Z, gₘₐₓ, kₚ, palatability)
end

end # module
