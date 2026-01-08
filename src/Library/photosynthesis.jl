"""Functions related to phytoplankton light uptake."""

module Photosynthesis

using Agate.Library.Nutrients: monod_limitation, liebig_minimum

export light_limitation_smith,
    light_limitation_geider,
    photosynthetic_growth_single_nutrient,
    photosynthetic_growth_single_nutrient_geider_light,
    photosynthetic_growth_two_nutrients_geider_light

"""
    light_limitation_smith(PAR, initial_slope, maximum_growth_0C)
    
Smith 1936 formulation of light limitation (also see Evans and Parslow, 1985).

!!! formulation
    α * PAR / √(μ₀² + α² * PAR²)

    where:
    - PAR = photosynthetic active radiation
    - α = initial photosynthetic slope
    - μ₀ = maximum growth rate at 0 °C    

# Arguments
- `PAR`: photosynthetic active radiation
- `initial_slope`: initial photosynthetic slope
- `maximum_growth_0C`: maximum growth rate at T = 0 °C
"""
@inline function light_limitation_smith(PAR, initial_slope, maximum_growth_0C)
    if initial_slope == zero(initial_slope)
        return zero(initial_slope)
    end

    return initial_slope * PAR / sqrt(maximum_growth_0C^2 + initial_slope^2 * PAR^2)
end

"""
    light_limitation_geider(PAR, photosynthetic_slope, maximum_growth_rate, 
        chlorophyll_to_carbon_ratio)

A light limitation function which depends on the cellular ratio of chlorophyll to carbon.
This formulation is based on equation (4) from Geider et al., 1998.

!!! formulation
    Pᶜₘₐₓ[1-exp((-αᶜʰˡθᶜPAR)/Pᶜₘₐₓ)]
    
    where:
    - Pᶜₘₐₓ = maximum_growth rate
    - αᶜʰˡ = photosynthetic slope
    - PAR = photosynthetic active radiation
    - θᶜ = chlorophyll to carbon ratio

# Arguments
- `PAR`: photosynthetic active radiation
- `maximum_growth_rate`: maximum growth rate before nutrient limitation
- `photosynthetic_slope`: initial photosynthetic slope
- `chlorophyll_to_carbon_ratio`: ratio between cellular chlorophyll and carbon
"""
@inline function light_limitation_geider(
    PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio
)
    if maximum_growth_rate == zero(maximum_growth_rate)
        return zero(maximum_growth_rate)
    end

    return maximum_growth_rate * (
        one(maximum_growth_rate) - exp(
            (-photosynthetic_slope * chlorophyll_to_carbon_ratio * PAR) /
            maximum_growth_rate,
        )
    )
end

"""
    photosynthetic_growth_single_nutrient(N, P, PAR, maximum_growth_0C, 
        nutrient_half_saturation, initial_slope)

Single nutrient monod smith photosynthetic growth (used, for example, in Kuhn 2015).

!!! formulation
    μ₀ * γᴿ * γᴾᴬᴿ * ``P``
    
    where:
    - μ₀ = maximum growth rate at 0 °C
    - γᴿ = `monod_limitation(R, kᵣ)`
    - R = nutrient concentration
    - kᵣ = nutrient half saturation
    - γᴾᴬᴿ = `light_limitation_smith(PAR, α, μ₀)`
    - PAR = photosynthetic active radiation
    - α = initial photosynthetic slope
    - P = plankton concentration

# Arguments
- `R`: nutrient concentration
- `P`: plankton concentration
- `PAR`: photosynthetic active radiation
- `maximum_growth_0C`: maximum growth rate at T = 0 °C
- `nutrient_half_saturation`: nutrient half saturation
- `initial_slope`: initial photosynthetic slope
"""
@inline function photosynthetic_growth_single_nutrient(
    R, P, PAR, maximum_growth_0C, nutrient_half_saturation, initial_slope
)
    return maximum_growth_0C *
           monod_limitation(R, nutrient_half_saturation) *
           light_limitation_smith(PAR, initial_slope, maximum_growth_0C) *
           P
end

"""
    photosynthetic_growth_single_nutrient_geider_light(R, P, PAR, maximum_growth_rate, 
        nutrient_half_saturation, photosynthetic_slope, chlorophyll_to_carbon_ratio)

Single nutrient geider photosynthetic growth.

!!! formulation
    γᴿ * γᴾᴬᴿ * P

    where:
    - γᴿ = `monod_limitation(R, kᵣ)`
    - R = nutrient concentration
    - kᵣ = nutrient half saturation
    - γᴾᴬᴿ = `light_limitation_geider(PAR, α, Pᶜₘₐₓ, θᶜ)`
    - PAR = photosynthetic active radiation
    - Pᶜₘₐₓ = maximum growth rate
    - θᶜ = chlorophyll to carbon ratio
    - α = initial photosynthetic slope
    - P = plankton concentration

# Arguments
- `N`: nutrient concentration
- `P`: phytoplankton concentration
- `PAR`: photosynthetic active radiation
- `maximum_growth_rate`: maximum growth rate before nutrient limitation
- `nutrient_half_saturation`: nutrient half saturation
- `photosynthetic_slope`: initial photosynthetic slope
- `chlorophyll_to_carbon_ratio`: ratio between cellular chlorophyll and carbon
"""
@inline function photosynthetic_growth_single_nutrient_geider_light(
    R,
    P,
    PAR,
    maximum_growth_rate,
    nutrient_half_saturation,
    photosynthetic_slope,
    chlorophyll_to_carbon_ratio,
)
    return monod_limitation(R, nutrient_half_saturation) *
           light_limitation_geider(
               PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio
           ) *
           P
end

"""
    photosynthetic_growth_two_nutrients_geider_light(
        DIN,
        PO4,
        P,
        PAR,
        maximum_growth_rate,
        half_saturation_DIN,
        half_saturation_PO4,
        photosynthetic_slope,
        chlorophyll_to_carbon_ratio,
    )

Two nutrient geider photosynthetic growth.


!!! formulation
    Pᶜₘ *  γᴾᴬᴿ * P

    where:
    - Pᶜₘ = `liebig_minimum([γᴺ, γᴾ])` * Pᶜₘₐₓ
    - γᴾᴬᴿ = `light_limitation_geider(PAR, α, Pᶜₘ, θᶜ)`
    - γᴺ = `monod_limitation(DIN, Kₙ)`
    - γᴾ = `monod_limitation(PO₄, Kₚ)`
    - DIN = dissolved inorganic nitrogen concentration
    - PO₄ = phosphate concentration
    - Kₙ = nitrogen half saturation constant
    - Kₚ = phosphate half saturation constant
    - Pᶜₘₐₓ = maximum growth rate
    - θᶜ = chlorophyll to carbon ratio
    - α = initial photosynthetic slope
    - P = plankton concentration

!!! info
    Unlike the MITgcm-DARWIN formulation this function does not include light inhibition

# Arguments
- `DIN`: dissolved inorganic nitrogen concentration
- `PO4`: phosphate concentration
- `P`: phytoplankton concentration
- `PAR`: photosynthetic active radiation
- `maximum_growth_rate`: maximum growth rate before nutrient limitation
- `half_saturation_DIN`: nitrogen half saturation
- `half_saturation_PO4`: phosphate half saturation
- `photosynthetic_slope`: initial photosynthetic slope
- `chlorophyll_to_carbon_ratio`: ratio between cellular chlorophyll and carbon
"""
@inline function photosynthetic_growth_two_nutrients_geider_light(
    R1,
    R2,
    P,
    PAR,
    maximum_growth_rate,
    nutrient_half_saturation_1,
    nutrient_half_saturation_2,
    photosynthetic_slope,
    chlorophyll_to_carbon_ratio,
)
    nutrient_factor = liebig_minimum(
        monod_limitation(R1, nutrient_half_saturation_1),
        monod_limitation(R2, nutrient_half_saturation_2),
    )

    return nutrient_factor *
           light_limitation_geider(
               PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio
           ) *
           P
end

end # module
