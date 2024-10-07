nutrient_limitation(N, kₙ) = N / (kₙ + N)

Q₁₀(T) = 1.88 ^ (T / 10) # T in °C

light_limitation(PAR, α, μ₀)= α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)

phytoplankton_metabolic_loss(P, lᵖⁿ) = lᵖⁿ  * P

zooplankton_metabolic_loss(Z, lᶻⁿ) = lᶻⁿ  * Z

remineralization(D, rᵈⁿ) = rᵈⁿ * D

phytoplankton_growth(N, P, PAR, μ₀, kₙ, α) = μ₀ * nutrient_limitation(N, kₙ) * light_limitation(PAR, α, μ₀) * P
phytoplankton_grazing_loss(P, Z, gₘₐₓ, kₚ) = gₘₐₓ * nutrient_limitation(P ^ 2, kₚ ^ 2) * Z
phytoplankton_mortality_loss(P, lᵖᵈ) = lᵖᵈ * P

zooplankton_grazing_gain(P, Z, β, gₘₐₓ, kₚ) = β * gₘₐₓ * nutrient_limitation(P ^ 2, kₚ ^ 2) * Z
zooplankton_mortality_loss(Z, lᶻᵈ) = lᶻᵈ * Z ^ 2
zooplankton_assimilation_loss(P, Z, β, gₘₐₓ, kₚ) = (1 - β) * gₘₐₓ * nutrient_limitation(P ^ 2, kₚ ^ 2) * Z

dummy_gmax_function(gₘₐₓ, a, b) = ((gₘₐₓ*a) + b-b)/a # this will always return gₘₐₓ
