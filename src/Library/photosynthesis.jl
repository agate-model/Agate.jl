"""Photosynthesis and light-uptake functors."""

module Photosynthesis

using ..Nutrients: MonodLimitation, LiebigMinimum

export SmithLightLimitation,
    GeiderLightLimitation,
    SingleNutrientGrowthSmith,
    SingleNutrientGrowthGeider,
    TwoNutrientGrowthGeider,
    smith_growth,
    geider_growth,
    geider_growth_two_nutrients

"""
    SmithLightLimitation(alpha, maximum_growth_0C)

Smith-style light limitation factor as a function of PAR.
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



"""
    smith_growth(maximum_growth_0C, nutrient_half_saturation, alpha, R, P, PAR)

Convenience wrapper for `SingleNutrientGrowthSmith(maximum_growth_0C, nutrient_half_saturation, alpha)(R, P, PAR)`.
"""
@inline smith_growth(maximum_growth_0C, nutrient_half_saturation, alpha, R, P, PAR) =
    SingleNutrientGrowthSmith(maximum_growth_0C, nutrient_half_saturation, alpha)(R, P, PAR)

"""
    geider_growth(maximum_growth_rate, nutrient_half_saturation, alpha, chlorophyll_to_carbon_ratio, R, P, PAR)

Convenience wrapper for `SingleNutrientGrowthGeider(maximum_growth_rate, nutrient_half_saturation, alpha, chlorophyll_to_carbon_ratio)(R, P, PAR)`.
"""
@inline geider_growth(maximum_growth_rate, nutrient_half_saturation, alpha, chlorophyll_to_carbon_ratio, R, P, PAR) =
    SingleNutrientGrowthGeider(maximum_growth_rate, nutrient_half_saturation, alpha, chlorophyll_to_carbon_ratio)(R, P, PAR)

"""
    geider_growth_two_nutrients(maximum_growth_rate, half_saturation_1, half_saturation_2, alpha, chlorophyll_to_carbon_ratio, R1, R2, P, PAR)

Convenience wrapper for `TwoNutrientGrowthGeider(maximum_growth_rate, half_saturation_1, half_saturation_2, alpha, chlorophyll_to_carbon_ratio)(R1, R2, P, PAR)`.
"""
@inline geider_growth_two_nutrients(
    maximum_growth_rate,
    half_saturation_1,
    half_saturation_2,
    alpha,
    chlorophyll_to_carbon_ratio,
    R1,
    R2,
    P,
    PAR,
) =
    TwoNutrientGrowthGeider(
        maximum_growth_rate,
        half_saturation_1,
        half_saturation_2,
        alpha,
        chlorophyll_to_carbon_ratio,
    )(R1, R2, P, PAR)

end # module
