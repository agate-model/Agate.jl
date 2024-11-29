"""
Functions related to zooplankton predation.
"""

module Predation

using ..Nutrients

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
holling_type_2(R::Real, k::Real) = R / (k + R)

"""
Estimates the loss rate of P (prey), to Z (predator).

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: grazing half saturation
"""
idealized_predation_loss(P, Z, gₘₐₓ, kₚ) = gₘₐₓ * monod_limitation(P^2, kₚ^2) * Z

"""
Estimates the gain rate of Z (predator) feeding on P (prey).

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: grazing half saturation
"""
idealized_predation_gain(P, Z, β, gₘₐₓ, kₚ) = β * gₘₐₓ * monod_limitation(P^2, kₚ^2) * Z

"""
Estimates the rate at which plankton predation gain is lost due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing half saturation
"""
function idealized_predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ)
    return (1 - β) * gₘₐₓ * monod_limitation(P^2, kₚ^2) * Z
end

end # module
