"""
Functions related to phytoplankton light uptake.
"""

module Photosynthesis

export γˡⁱᵍʰᵗ, smith_light_limitation, idealized_photosynthetic_growth

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
    if alpha == 0
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
- `α`: initial photosynthetic slope
- `μ₀`: maximum growth rate at T = 0 °C
- `kₙ`: nutrient half saturation
"""
function idealized_photosynthetic_growth(N, P, PAR, μ₀, kₙ, α)
    return μ₀ * monod_limitation(N, kₙ) * smith_light_limitation(PAR, α, μ₀) * P
end

end # module
