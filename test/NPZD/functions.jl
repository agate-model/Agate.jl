monod_limitation(R, k) = R / (k + R)
holling_type_2(R, k) = R / (k + R)

smith_light_limitation(PAR, α, μ₀) = α * PAR / sqrt(μ₀^2 + α^2 * PAR^2)

remineralization_idealized(D, r) = r * D

function photosynthetic_growth_idealized(N, P, PAR, μ₀, kₙ, α)
    return μ₀ * monod_limitation(N, kₙ) * smith_light_limitation(PAR, α, μ₀) * P
end

predation_loss_idealized(P, Z, gₘₐₓ, kₚ) = gₘₐₓ * holling_type_2(P^2, kₚ^2) * Z

predation_gain_idealized(P, Z, β, gₘₐₓ, kₚ) = β * gₘₐₓ * holling_type_2(P^2, kₚ^2) * Z

function predation_assimilation_loss_idealized(P, Z, β, gₘₐₓ, kₚ)
    return (1 - β) * gₘₐₓ * holling_type_2(P^2, kₚ^2) * Z
end

linear_loss(P, l) = l * P

quadratic_loss(P, l) = l * P^2
