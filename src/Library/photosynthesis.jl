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

Smith (1936) light limitation formulation.

Returns zero when `initial_slope` is zero.
"""
@inline function light_limitation_smith(PAR, initial_slope, maximum_growth_0C)
    if initial_slope == zero(initial_slope)
        return zero(initial_slope)
    end

    return initial_slope * PAR / sqrt(maximum_growth_0C^2 + initial_slope^2 * PAR^2)
end

"""
    light_limitation_geider(PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio)

Geider-style light limitation that depends on chlorophyll-to-carbon ratio.

Returns zero when `maximum_growth_rate` is zero.
"""
@inline function light_limitation_geider(
    PAR,
    photosynthetic_slope,
    maximum_growth_rate,
    chlorophyll_to_carbon_ratio,
)
    if maximum_growth_rate == zero(maximum_growth_rate)
        return zero(maximum_growth_rate)
    end

    return maximum_growth_rate * (
        one(maximum_growth_rate) -
        exp((-photosynthetic_slope * chlorophyll_to_carbon_ratio * PAR) / maximum_growth_rate)
    )
end

"""
    photosynthetic_growth_single_nutrient(R, P, PAR, maximum_growth_0C, nutrient_half_saturation, initial_slope)

Monod nutrient limitation combined with Smith light limitation.
"""
@inline function photosynthetic_growth_single_nutrient(
    R,
    P,
    PAR,
    maximum_growth_0C,
    nutrient_half_saturation,
    initial_slope,
)
    return maximum_growth_0C *
           monod_limitation(R, nutrient_half_saturation) *
           light_limitation_smith(PAR, initial_slope, maximum_growth_0C) *
           P
end

"""
    photosynthetic_growth_single_nutrient_geider_light(
        R,
        P,
        PAR,
        maximum_growth_rate,
        nutrient_half_saturation,
        photosynthetic_slope,
        chlorophyll_to_carbon_ratio,
    )

Monod nutrient limitation combined with Geider-style light limitation.
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
           light_limitation_geider(PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio) *
           P
end

"""
    photosynthetic_growth_two_nutrients_geider_light(
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

Two-nutrient growth using Liebig's law of the minimum, combined with Geider-style light limitation.
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
           light_limitation_geider(PAR, photosynthetic_slope, maximum_growth_rate, chlorophyll_to_carbon_ratio) *
           P
end

end # module
