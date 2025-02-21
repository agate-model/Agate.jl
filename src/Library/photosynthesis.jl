"""
Functions related to phytoplankton light uptake.
"""

module Photosynthesis

using Agate.Library.Nutrients

export γˡⁱᵍʰᵗ,
    light_limitation_smith,
    light_limitation_geider,
    photosynthetic_growth_single_nutrient,
    photosynthetic_growth_single_nutrient_geider_light,
    net_photosynthetic_growth_single_nutrient,
    net_photosynthetic_growth_single_nutrient_geider_light
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
function light_limitation_smith(PAR, α, μ₀)
    # here to avoid division by 0 when α and μ₀ are both 0
    if α == 0
        return 0.0
    end
    return α * PAR / sqrt(μ₀^2 + α^2 * PAR^2)
end

"""
    Pᶜₘₐₓ[1-exp((-αᶜʰˡθᶜE₀)/Pᶜₘₐₓ)]

A light limitation function which depends on the cellular ratio of chlorophyll to carbon.
This formulation is based on equation (4) from Geider et al., 1998.

# Arguments
- `PAR`: photosynthetic active radiation (E₀)
- `maximum_growth_rate`: maximum growth rate before nutrient limitation (Pᶜₘₐₓ)
- `photosynthetic_slope`: initial photosynthetic slope (αᶜʰˡ)
- `chlorophyll_to_carbon_ratio`: ratio between cellular chlorophyll and carbon (θᶜ)
"""
function light_limitation_geider(
    PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio
)
    if maximum_growth_rate == 0
        return 0.0
    end
    return maximum_growth_rate * (
        1 - exp(
            (-photosynthetic_slope * chlorophyll_to_carbon_ratio * PAR) /
            maximum_growth_rate,
        )
    )
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
function photosynthetic_growth_single_nutrient(N, P, PAR, μ₀, kₙ, α)
    return μ₀ * monod_limitation(N, kₙ) * light_limitation_smith(PAR, α, μ₀) * P
end

"""
Single nutrient geider photosynthetic growth.

# Arguments
- `N`: nutrient concentration
- `P`: phytoplankton concentration
- `PAR`: photosynthetic active radiation
- `maximum_growth_rate`: maximum growth rate before nutrient limitation (Pᶜₘₐₓ)
- `kₙ`: nutrient half saturation
- `photosynthetic_slope`: initial photosynthetic slope (αᶜʰˡ)
- `chlorophyll_to_carbon_ratio`: ratio between cellular chlorophyll and carbon (θᶜ)
"""
function photosynthetic_growth_single_nutrient_geider_light(
    N, P, PAR, maximum_growth_rate, kₙ, photosynthetic_slope, chlorophyll_to_carbon_ratio
)
    return monod_limitation(N, kₙ) *
           light_limitation_geider(
               PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio
           ) *
           P
end

"""
Net photosynthetic growth of all plankton.

# Arguments
- `N`: Nutrient
- `P`: NamedArray which includes all plankton concentration values
- `PAR`: PAR
- `maximum_growth_rate`: NamedArray of all plankton maximum growth rates
- `nutrient_half_saturation`: NamedArray of all plankton nutrient half saturation constants
- `photosynthetic_slope`: initial photosynthetic slope (αᶜʰˡ)
- `chlorophyll_to_carbon_ratio`: ratio between cellular chlorophyll and carbon (θᶜ)
"""
function net_photosynthetic_growth_single_nutrient_geider_light(
    N,
    P,
    PAR,
    maximum_growth_rate,
    nutrient_half_saturation,
    photosynthetic_slope,
    chlorophyll_to_carbon_ratio,
)
    sum(photosynthetic_growth_single_nutrient_geider_light.(N, P, PAR, maximum_growth_rate, nutrient_half_saturation, photosynthetic_slope, chlorophyll_to_carbon_ratio))
end

"""
Net photosynthetic growth of all plankton assuming geider light limitation.

# Arguments
- `N`: Nutrient
- `P`: NamedArray which includes all plankton concentration values
- `PAR`: PAR
- `maximum_growth_rate`: NamedArray of all plankton maximum growth rates
- `nutrient_half_saturation`: NamedArray of all plankton nutrient half saturation constants
- `alpha`: initial photosynthetic slope
"""
function net_photosynthetic_growth_single_nutrient(
    N,
    P,
    PAR,
    maximum_growth_rate,
    nutrient_half_saturation,
    alpha,
)
    sum(photosynthetic_growth_single_nutrient.(N, P, PAR, maximum_growth_rate, nutrient_half_saturation, alpha))
end

end # module
