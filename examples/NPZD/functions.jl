menden_limitation(R, k) = R / (k + R)

Q₁₀(T) = 1.88^(T / 10) # T in °C

smith_light_limitation(PAR, α, μ₀) = α * PAR / sqrt(μ₀^2 + α^2 * PAR^2)

remineralization(D, r) = r * D

function photosynthetic_growth(N, P, PAR, μ₀, kₙ, α)
    return μ₀ * menden_limitation(N, kₙ) * smith_light_limitation(PAR, α, μ₀) * P
end

predation_loss(P, Z, gₘₐₓ, kₚ) = gₘₐₓ * menden_limitation(P^2, kₚ^2) * Z

predation_gain(P, Z, β, gₘₐₓ, kₚ) = β * gₘₐₓ * menden_limitation(P^2, kₚ^2) * Z

function predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ)
    return (1 - β) * gₘₐₓ * menden_limitation(P^2, kₚ^2) * Z
end

linear_loss(P, l) = l * P

quadratic_loss(P, l) = l * P^2

function summed_linear_loss(P, l)
    return sum([
        linear_loss(P[plankton_index], l[plankton_index]) for plankton_index in eachindex(P)
    ])
end
