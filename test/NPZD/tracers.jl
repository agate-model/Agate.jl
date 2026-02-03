using Agate
using Agate.Functors: CompiledEquation, Requirements
using Agate.Constructor: define_tracer_functions
using Oceananigans.Units

include(joinpath(@__DIR__, "functions.jl"))

# Test NPZD model used to validate that Agate-generated tracer functions integrate
# consistently with OceanBioME's BoxModel infrastructure.

parameters = (;
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

fN =
    (bgc, x, y, z, t, N, D, P, Z, PAR) -> begin
        p = bgc.parameters
        linear_loss(P, p.lᵖⁿ) +
        linear_loss(Z, p.lᶻⁿ) +
        remineralization_idealized(D, p.rᵈⁿ) -
        photosynthetic_growth_idealized(N, P, PAR, p.μ₀, p.kₙ, p.α)
    end

fD =
    (bgc, x, y, z, t, N, D, P, Z, PAR) -> begin
        p = bgc.parameters
        linear_loss(P, p.lᵖᵈ) +
        predation_assimilation_loss_idealized(P, Z, p.β, p.gₘₐₓ, p.kₚ) +
        quadratic_loss(Z, p.lᶻᵈ) - remineralization_idealized(D, p.rᵈⁿ)
    end

fP =
    (bgc, x, y, z, t, N, D, P, Z, PAR) -> begin
        p = bgc.parameters
        photosynthetic_growth_idealized(N, P, PAR, p.μ₀, p.kₙ, p.α) -
        predation_loss_idealized(P, Z, p.gₘₐₓ, p.kₚ) - linear_loss(P, p.lᵖⁿ) -
        linear_loss(P, p.lᵖᵈ)
    end

fZ =
    (bgc, x, y, z, t, N, D, P, Z, PAR) -> begin
        p = bgc.parameters
        predation_gain_idealized(P, Z, p.β, p.gₘₐₓ, p.kₚ) - linear_loss(Z, p.lᶻⁿ) -
        quadratic_loss(Z, p.lᶻᵈ)
    end

tracers = (
    N=CompiledEquation(fN, Requirements(; scalars=(:lᵖⁿ, :lᶻⁿ, :rᵈⁿ, :μ₀, :kₙ, :α))),
    D=CompiledEquation(fD, Requirements(; scalars=(:lᵖᵈ, :β, :gₘₐₓ, :kₚ, :lᶻᵈ, :rᵈⁿ))),
    P=CompiledEquation(fP, Requirements(; scalars=(:μ₀, :kₙ, :α, :gₘₐₓ, :kₚ, :lᵖⁿ, :lᵖᵈ))),
    Z=CompiledEquation(fZ, Requirements(; scalars=(:β, :gₘₐₓ, :kₚ, :lᶻⁿ, :lᶻᵈ))),
)

AgateNPZD = define_tracer_functions(parameters, tracers)
