"""
Functions related to phytoplankton light uptake.
"""

module Photosynthesis

using Agate.Library.Nutrients

export γˡⁱᵍʰᵗ,
    smith_light_limitation,
    idealized_photosynthetic_growth,
    net_idealized_photosynthetic_growth

"""
    γˡⁱᵍʰᵗ = (1 - ℯ^(kˢᵃᵗ*I)) * ℯ^kⁱⁿʰ * nˡⁱᵍʰᵗ

Light limitation for plankton (Default MITgcm-DARWIN formulation).

# Arguments
- `I`: irradiance
- `kˢᵃᵗ`:  half saturation constant of light saturation
- `kⁱⁿʰ`: half saturation constant of light inhibition
- `nˡⁱᵍʰᵗ`: light penalty term
"""
function γˡⁱᵍʰᵗ(I, kˢᵃᵗ, kⁱⁿʰ, nˡⁱᵍʰᵗ)
    return (1 - ℯ^(kˢᵃᵗ * I)) * ℯ^kⁱⁿʰ * nˡⁱᵍʰᵗ
end

"
    α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)

Smith 1936 formulation of light limitation (also see Evans and Parslow, 1985).

# Arguments
- `PAR`: photosynthetic active radiation
- `α`: initial photosynthetic slope
- `μ₀`: maximum growth rate at T = 0 °C (this seems weird?, from Kuhn 2015)
"
function smith_light_limitation(PAR, α, μ₀)
    # here to avoid division by 0 when α and μ₀ are both 0
    if α == 0
        return 0.0
    end
    return α * PAR / sqrt(μ₀^2 + α^2 * PAR^2)
end

"""
Single nutrient monod smith photosynthetic growth (used, for example, in Kuhn 2015).

# Arguments
- `N`: nutrient concentration
- `P`: phytoplankton concentration
- `PAR`: photosynthetic active radiation
- `μ₀`: maximum growth rate at T = 0 °C
- `kₙ`: nutrient half saturation
- `α`: initial photosynthetic slope
"""
function idealized_photosynthetic_growth(N, P, PAR, μ₀, kₙ, α)
    return μ₀ * monod_limitation(N, kₙ) * smith_light_limitation(PAR, α, μ₀) * P
end

"""
Net photosynthetic growth of all plankton.

# Arguments
- `N`: Nitrogen
- `P`: NamedArray which includes all plankton
- `PAR`: PAR
- `maximum_growth_rate`: NamedArray of all plankton maximum growth rates
- `nitrogen_half_saturation`: NamedArray of all plankton nitrogen half saturation constants
"""
function net_idealized_photosynthetic_growth(
    N, P, PAR, maximum_growth_rate, nitrogen_half_saturation, alpha
)
    return sum([
        # only phytoplankton have maximum_growth_rate, nitrogen_half_saturation and alpha
        # --> get names from either of those arrays
        idealized_photosynthetic_growth(
            N,
            P[name],
            PAR,
            maximum_growth_rate[name],
            nitrogen_half_saturation[name],
            alpha[name],
        ) for name in names(maximum_growth_rate)[1]
    ],)
end

end # module
