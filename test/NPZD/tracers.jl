using Agate
using Oceananigans.Units

"""Runtime parameters for the test NPZD model."""
struct NPZDParameters{FT<:AbstractFloat}
    μ₀::FT
    kₙ::FT
    lᵖⁿ::FT
    lᶻⁿ::FT
    lᵖᵈ::FT
    gₘₐₓ::FT
    kₚ::FT
    β::FT
    lᶻᵈ::FT
    rᵈⁿ::FT
    α::FT
end

parameters = NPZDParameters{Float64}(
    0.6989 / day,
    2.3868,
    0.066 / day,
    0.0102 / day,
    0.0101 / day,
    2.1522 / day,
    0.5573,
    0.9116,
    0.3395 / day,
    0.1213 / day,
    0.1953 / day,
)

tracers = (
    N=:(
        linear_loss(P, lᵖⁿ) + linear_loss(Z, lᶻⁿ) + remineralization_idealized(D, rᵈⁿ) -
        photosynthetic_growth_idealized(N, P, PAR, μ₀, kₙ, α)
    ),
    D=:(
        linear_loss(P, lᵖᵈ) +
        predation_assimilation_loss_idealized(P, Z, β, gₘₐₓ, kₚ) +
        quadratic_loss(Z, lᶻᵈ) - remineralization_idealized(D, rᵈⁿ)
    ),
    P=:(
        photosynthetic_growth_idealized(N, P, PAR, μ₀, kₙ, α) -
        predation_loss_idealized(P, Z, gₘₐₓ, kₚ) - linear_loss(P, lᵖⁿ) - linear_loss(P, lᵖᵈ)
    ),
    Z=:(
        predation_gain_idealized(P, Z, β, gₘₐₓ, kₚ) - linear_loss(Z, lᶻⁿ) -
        quadratic_loss(Z, lᶻᵈ)
    ),
)

NPZD = define_tracer_functions(
    parameters,
    tracers;
    helper_functions=joinpath(@__DIR__, "functions.jl"),
)
