"""
Functions related to phytoplankton light uptake
"""

module Photosynthesis

export γⱼˡⁱᵍʰᵗ, smith_light_limitation, idealized_photosynthetic_growth

# TODO: get rid of j subscript in the function args
"
    γⱼˡⁱᵍʰᵗ = (1 - ℯ^(kⱼˢᵃᵗ*I)) * ℯ^kⱼⁱⁿʰ * nⱼˡⁱᵍʰᵗ

Light limitation for plankton j (Default MITgcm-DARWIN formulation).

# Arguments
- `I`: irradiance,
- `kⱼˢᵃᵗ`:  half saturation constant of light saturation of plankton j,
- `kⱼⁱⁿʰ`: half saturation constant of light inhibition of plankton j,
- `nⱼˡⁱᵍʰᵗ`: light penalty term of plankton j
"
function γⱼˡⁱᵍʰᵗ(I, kⱼˢᵃᵗ, kⱼⁱⁿʰ, nⱼˡⁱᵍʰᵗ)
    return (1 - ℯ^(kⱼˢᵃᵗ * I)) * ℯ^kⱼⁱⁿʰ * nⱼˡⁱᵍʰᵗ
end

"
    α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)

Smith 1936 formulation of light limitation (also see Evans and Parslow, 1985).

# Arguments
- `PAR`: Photosynthetic Active Radiation
- `α`: Initial photosynthetic slope
- `μ₀`: Maximum growth rate at T = 0 °C (this seems weird?, from Kuhn 2015)
"
function smith_light_limitation(PAR, α, μ₀)
    # here to avoid division by 0 when α and μ₀ are both 0
    if alpha == 0
        return 0.0
    end
    return α * PAR / sqrt(μ₀^2 + α^2 * PAR^2)
end

"""
Single nutrient monod smith photosynthetic growth (used, for example, in Kuhn 2015)

# Arguments
- `N`: nutrients concentration ?
- `P`: phytoplankton concentration ?
- `PAR`: Photosynthetic Active Radiation
- `α`: Initial photosynthetic slope
- `μ₀`: Maximum growth rate at T = 0 °C
- `kₙ`: ??
"""
function idealized_photosynthetic_growth(N, P, PAR, μ₀, kₙ, α)
    return μ₀ * menden_limitation(N, kₙ) * smith_light_limitation(PAR, α, μ₀) * P
end

end # module
