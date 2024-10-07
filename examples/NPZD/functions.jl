menden_limitation(R, k) = R / (k + R)

Q₁₀(T) = 1.88 ^ (T / 10) # T in °C

smith_light_limitation(PAR, α, μ₀)= α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)

remineralization(D, r) = r * D

photosynthetic_growth(N, P, PAR, μ₀, kₙ, α) = μ₀ * menden_limitation(N, kₙ) * smith_light_limitation(PAR, α, μ₀) * P

predation_loss(P, Z, gₘₐₓ, kₚ) = gₘₐₓ * menden_limitation(P ^ 2, kₚ ^ 2) * Z

predation_gain(P, Z, β, gₘₐₓ, kₚ) = β * gₘₐₓ * menden_limitation(P ^ 2, kₚ ^ 2) * Z

predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ) = (1 - β) * gₘₐₓ * menden_limitation(P ^ 2, kₚ ^ 2) * Z

linear_loss(P, l) = l  * P

quadratic_loss(P, l) = l  * P ^ 2
