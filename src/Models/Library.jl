nutrient_limitation(N; kₙ=kₙ) = N / (kₙ + N)

Q₁₀(T) = 1.88 ^ (T / 10) # T in °C

light_limitation(PAR; α=α, μ₀=μ₀)= α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)

phytoplankton_metabolic_loss(P; lᵖⁿ=lᵖⁿ) = lᵖⁿ  * P

zooplankton_metabolic_loss(Z; lᶻⁿ=lᶻⁿ) = lᶻⁿ  * Z

remineralization(D; rᵈⁿ=rᵈⁿ) = rᵈⁿ * D

phytoplankton_growth(N, P, PAR; μ₀=μ₀, kₙ=kₙ, α=α) = μ₀ * nutrient_limitation(N, kₙ=kₙ) * light_limitation(PAR, α=α, μ₀=μ₀) * P
phytoplankton_grazing_loss(P, Z; gₘₐₓ=gₘₐₓ, kₚ=kₚ) = gₘₐₓ * nutrient_limitation(P ^ 2, kₙ=kₚ ^ 2) * Z
phytoplankton_mortality_loss(P; lᵖᵈ=lᵖᵈ) = lᵖᵈ * P

zooplankton_grazing_gain(P, Z; β=β, gₘₐₓ=gₘₐₓ, kₚ=kₚ) = β * gₘₐₓ * nutrient_limitation(P ^ 2, kₙ=kₚ ^ 2) * Z
zooplankton_mortality_loss(Z; lᶻᵈ=lᶻᵈ) = lᶻᵈ * Z ^ 2
zooplankton_assimilation_loss(P, Z; β=β, gₘₐₓ=gₘₐₓ, kₚ=kₚ) = (1 - β) * gₘₐₓ * nutrient_limitation(P ^ 2, kₙ=kₚ ^ 2) * Z
