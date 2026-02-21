"""Photosynthesis and light-uptake functors."""

module Photosynthesis

using ..Nutrients: MonodLimitation, LiebigMinimum

export smith_light_limitation,
    geider_light_limitation,
    smith_single_nutrient_growth,
    geider_single_nutrient_growth,
    geider_two_nutrient_growth

"""
    SmithLightLimitation(alpha, maximum_growth_0C)

Smith (1936) light limitation factor as a function of PAR.

!!! formulation
    α * PAR / √(μ₀² + α² * PAR²)

    where:
    - PAR = photosynthetically active radiation
    - α = initial photosynthetic slope
    - μ₀ = maximum growth rate at 0 °C
"""
struct SmithLightLimitation{T}
    alpha::T
    maximum_growth_0C::T
end

@inline function (f::SmithLightLimitation)(PAR)
    α = f.alpha
    μ₀ = f.maximum_growth_0C
    if α == zero(α) || μ₀ == zero(μ₀)
        return zero(α)
    end
    return α * PAR / sqrt(μ₀ * μ₀ + α * α * PAR * PAR)
end

"""
    GeiderLightLimitation(alpha, maximum_growth_rate, chlorophyll_to_carbon_ratio)

Geider-style light limitation factor as a function of PAR.

This formulation follows Geider et al. (1998), equation (4).

!!! formulation
    Pᶜₘₐₓ[1-exp((-αᶜʰˡθᶜPAR)/Pᶜₘₐₓ)]

    where:
    - Pᶜₘₐₓ = maximum growth rate
    - αᶜʰˡ = photosynthetic slope
    - PAR = photosynthetically active radiation
    - θᶜ = chlorophyll to carbon ratio
"""
struct GeiderLightLimitation{T}
    alpha::T
    maximum_growth_rate::T
    chlorophyll_to_carbon_ratio::T
end

@inline function (f::GeiderLightLimitation)(PAR)
    α = f.alpha
    Pᶜₘₐₓ = f.maximum_growth_rate
    θᶜ = f.chlorophyll_to_carbon_ratio
    if Pᶜₘₐₓ == zero(Pᶜₘₐₓ)
        return zero(Pᶜₘₐₓ)
    end
    return Pᶜₘₐₓ * (one(Pᶜₘₐₓ) - exp((-α * θᶜ * PAR) / Pᶜₘₐₓ))
end

"""
    SingleNutrientGrowthSmith(maximum_growth_0C, nutrient_half_saturation, alpha)

Photosynthetic growth with single-nutrient Monod limitation and Smith light limitation.

!!! formulation
    μ₀ * γᴿ * γᴾᴬᴿ * P

    where:
    - μ₀ = maximum growth rate at 0 °C
    - γᴿ = `MonodLimitation(kᵣ)(R)`
    - R = nutrient concentration
    - kᵣ = nutrient half-saturation
    - γᴾᴬᴿ = `SmithLightLimitation(α, μ₀)(PAR)`
    - PAR = photosynthetically active radiation
    - α = initial photosynthetic slope
    - P = plankton concentration
"""
struct SingleNutrientGrowthSmith{T}
    maximum_growth_0C::T
    nutrient_half_saturation::T
    alpha::T
end

@inline function (g::SingleNutrientGrowthSmith)(R, P, PAR)
    μ₀ = g.maximum_growth_0C
    return μ₀ *
           MonodLimitation(g.nutrient_half_saturation)(R) *
           SmithLightLimitation(g.alpha, μ₀)(PAR) *
           P
end

"""
    SingleNutrientGrowthGeider(maximum_growth_rate, nutrient_half_saturation, alpha, chlorophyll_to_carbon_ratio)

Photosynthetic growth with single-nutrient Monod limitation and Geider light limitation.

!!! formulation
    γᴿ * γᴾᴬᴿ * P

    where:
    - γᴿ = `MonodLimitation(kᵣ)(R)`
    - R = nutrient concentration
    - kᵣ = nutrient half-saturation
    - γᴾᴬᴿ = `GeiderLightLimitation(α, Pᶜₘₐₓ, θᶜ)(PAR)`
    - PAR = photosynthetically active radiation
    - Pᶜₘₐₓ = maximum growth rate
    - θᶜ = chlorophyll to carbon ratio
    - α = photosynthetic slope
    - P = plankton concentration
"""
struct SingleNutrientGrowthGeider{T}
    maximum_growth_rate::T
    nutrient_half_saturation::T
    alpha::T
    chlorophyll_to_carbon_ratio::T
end

@inline function (g::SingleNutrientGrowthGeider)(R, P, PAR)
    return MonodLimitation(g.nutrient_half_saturation)(R) *
           GeiderLightLimitation(
               g.alpha, g.maximum_growth_rate, g.chlorophyll_to_carbon_ratio
           )(
               PAR
           ) *
           P
end

"""
    TwoNutrientGrowthGeider(maximum_growth_rate, half_saturation_1, half_saturation_2, alpha, chlorophyll_to_carbon_ratio)

Photosynthetic growth with two-nutrient Liebig limitation and Geider light limitation.

!!! formulation
    Pᶜₘ * γᴾᴬᴿ * P

    where:
    - Pᶜₘ = `LiebigMinimum()(γ¹, γ²)` * Pᶜₘₐₓ
    - γᴾᴬᴿ = `GeiderLightLimitation(α, Pᶜₘ, θᶜ)(PAR)`
    - γ¹ = `MonodLimitation(K₁)(R₁)`
    - γ² = `MonodLimitation(K₂)(R₂)`
    - R₁, R₂ = nutrient concentrations
    - K₁, K₂ = nutrient half-saturation constants
    - Pᶜₘₐₓ = maximum growth rate
    - θᶜ = chlorophyll to carbon ratio
    - α = photosynthetic slope
    - P = plankton concentration

!!! info
    Unlike the MITgcm-DARWIN formulation, this growth function does not include light inhibition.
"""
struct TwoNutrientGrowthGeider{T}
    maximum_growth_rate::T
    half_saturation_1::T
    half_saturation_2::T
    alpha::T
    chlorophyll_to_carbon_ratio::T
end

@inline function (g::TwoNutrientGrowthGeider)(R1, R2, P, PAR)
    γ = LiebigMinimum()(
        MonodLimitation(g.half_saturation_1)(R1), MonodLimitation(g.half_saturation_2)(R2)
    )
    return γ *
           GeiderLightLimitation(
               g.alpha, g.maximum_growth_rate, g.chlorophyll_to_carbon_ratio
           )(
               PAR
           ) *
           P
end

# -----------------------------------------------------------------------------
# Explicit function aliases (preferred developer UX).
# -----------------------------------------------------------------------------

"""
    smith_light_limitation(PAR, alpha, maximum_growth_0C)

Smith (1936) formulation of light limitation.

!!! formulation
    α * PAR / √(μ₀² + α² * PAR²)

# Arguments
- `PAR`: photosynthetically active radiation
- `alpha`: initial photosynthetic slope α
- `maximum_growth_0C`: maximum growth rate μ₀ at T = 0 °C
"""
@inline smith_light_limitation(PAR, alpha, maximum_growth_0C) =
    SmithLightLimitation(alpha, maximum_growth_0C)(PAR)

"""
    geider_light_limitation(PAR, alpha, maximum_growth_rate, chlorophyll_to_carbon_ratio)

Geider-style light limitation.

!!! formulation
    Pᶜₘₐₓ[1-exp((-αᶜʰˡθᶜPAR)/Pᶜₘₐₓ)]

# Arguments
- `PAR`: photosynthetically active radiation
- `alpha`: photosynthetic slope αᶜʰˡ
- `maximum_growth_rate`: maximum growth rate Pᶜₘₐₓ
- `chlorophyll_to_carbon_ratio`: chlorophyll-to-carbon ratio θᶜ
"""
@inline geider_light_limitation(
    PAR, alpha, maximum_growth_rate, chlorophyll_to_carbon_ratio
) = GeiderLightLimitation(alpha, maximum_growth_rate, chlorophyll_to_carbon_ratio)(PAR)

"""
    smith_single_nutrient_growth(R, P, PAR, maximum_growth_0C, nutrient_half_saturation, alpha)

Single-nutrient photosynthetic growth with Monod nutrient limitation and Smith light limitation.

!!! formulation
    μ₀ * γᴿ * γᴾᴬᴿ * P

# Arguments
- `R`: nutrient concentration
- `P`: plankton concentration
- `PAR`: photosynthetically active radiation
- `maximum_growth_0C`: maximum growth rate μ₀ at T = 0 °C
- `nutrient_half_saturation`: nutrient half-saturation kᵣ
- `alpha`: initial photosynthetic slope α
"""
@inline smith_single_nutrient_growth(
    R, P, PAR, maximum_growth_0C, nutrient_half_saturation, alpha
) = SingleNutrientGrowthSmith(maximum_growth_0C, nutrient_half_saturation, alpha)(R, P, PAR)

"""
    geider_single_nutrient_growth(R, P, PAR, maximum_growth_rate, nutrient_half_saturation, alpha, chlorophyll_to_carbon_ratio)

Single-nutrient photosynthetic growth with Monod nutrient limitation and Geider light limitation.

!!! formulation
    γᴿ * γᴾᴬᴿ * P

# Arguments
- `R`: nutrient concentration
- `P`: plankton concentration
- `PAR`: photosynthetically active radiation
- `maximum_growth_rate`: maximum growth rate Pᶜₘₐₓ
- `nutrient_half_saturation`: nutrient half-saturation kᵣ
- `alpha`: photosynthetic slope α
- `chlorophyll_to_carbon_ratio`: chlorophyll-to-carbon ratio θᶜ
"""
@inline geider_single_nutrient_growth(
    R,
    P,
    PAR,
    maximum_growth_rate,
    nutrient_half_saturation,
    alpha,
    chlorophyll_to_carbon_ratio,
) = SingleNutrientGrowthGeider(
    maximum_growth_rate, nutrient_half_saturation, alpha, chlorophyll_to_carbon_ratio
)(
    R, P, PAR
)

"""
    geider_two_nutrient_growth(R1, R2, P, PAR, maximum_growth_rate, half_saturation_1, half_saturation_2, alpha, chlorophyll_to_carbon_ratio)

Two-nutrient photosynthetic growth with Liebig limitation and Geider light limitation.

!!! formulation
    Pᶜₘ * γᴾᴬᴿ * P

# Arguments
- `R1`: first nutrient concentration
- `R2`: second nutrient concentration
- `P`: plankton concentration
- `PAR`: photosynthetically active radiation
- `maximum_growth_rate`: maximum growth rate Pᶜₘₐₓ
- `half_saturation_1`: half-saturation K₁
- `half_saturation_2`: half-saturation K₂
- `alpha`: photosynthetic slope α
- `chlorophyll_to_carbon_ratio`: chlorophyll-to-carbon ratio θᶜ

!!! info
    Unlike the MITgcm-DARWIN formulation, this growth function does not include light inhibition.
"""
@inline geider_two_nutrient_growth(
    R1,
    R2,
    P,
    PAR,
    maximum_growth_rate,
    half_saturation_1,
    half_saturation_2,
    alpha,
    chlorophyll_to_carbon_ratio,
) = TwoNutrientGrowthGeider(
    maximum_growth_rate,
    half_saturation_1,
    half_saturation_2,
    alpha,
    chlorophyll_to_carbon_ratio,
)(
    R1, R2, P, PAR
)

end # module
