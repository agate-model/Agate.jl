
using Agate
using Oceananigans.Units

parameters = (
    μ₀=0.6989 / day,
    kₙ=2.3868,
    lᵖⁿ=0.066 / day,
    lᶻⁿ=0.0102 / day,
    lᵖᵈ=0.0101 / day,
    gₘₐₓ=2.1522 / day,
    kₚ=0.5573,
    β=0.9116,
    lᶻᵈ=0.3395 / day,
    rᵈⁿ=0.1213 / day,
    α=0.1953 / day,
)

tracers = Dict(
    "N" => :(
        linear_loss(P, lᵖⁿ) + linear_loss(Z, lᶻⁿ) + remineralization_idealized(D, rᵈⁿ) -
        photosynthetic_growth_idealized(N, P, PAR, μ₀, kₙ, α)
    ),
    "D" => :(
        linear_loss(P, lᵖᵈ) +
        predation_assimilation_loss_idealized(P, Z, β, gₘₐₓ, kₚ) +
        quadratic_loss(Z, lᶻᵈ) - remineralization_idealized(D, rᵈⁿ)
    ),
    "P" => :(
        photosynthetic_growth_idealized(N, P, PAR, μ₀, kₙ, α) -
        predation_loss_idealized(P, Z, gₘₐₓ, kₚ) - linear_loss(P, lᵖⁿ) -
        linear_loss(P, lᵖᵈ)
    ),
    "Z" => :(
        predation_gain_idealized(P, Z, β, gₘₐₓ, kₚ) - linear_loss(Z, lᶻⁿ) -
        quadratic_loss(Z, lᶻᵈ)
    ),
)

NPZD = define_tracer_functions(parameters, tracers; helper_functions="functions.jl")
