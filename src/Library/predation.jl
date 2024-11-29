module Predation

using ..Nutrients

export holling_type_2,
    idealized_predation_loss,
    idealized_predation_gain,
    idealized_predation_assimilation_loss

"""
Zookplankton growth.
"""
holling_type_2(R::Real, k::Real) = R / (k + R)

idealized_predation_loss(P, Z, gₘₐₓ, kₚ) = gₘₐₓ * menden_limitation(P^2, kₚ^2) * Z

idealized_predation_gain(P, Z, β, gₘₐₓ, kₚ) = β * gₘₐₓ * menden_limitation(P^2, kₚ^2) * Z

function idealized_predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ)
    return (1 - β) * gₘₐₓ * menden_limitation(P^2, kₚ^2) * Z
end

end # module
