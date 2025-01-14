module Tracers

using NamedArrays

export typical_detritus,
    typical_nutrients, simplified_phytoplankton_growth, simplified_zooplankton_growth

"""
Build expression for a single nutrient function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function typical_nutrients(plankton_array)
    return :(
        net_linear_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            linear_mortality,
            mortality_export_fraction,
        ) +
        net_quadratic_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            quadratic_mortality,
            mortality_export_fraction,
        ) +
        idealized_remineralization(D, detritus_remineralization) -
        net_idealized_photosynthetic_growth(
            N,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            PAR,
            maximum_growth_rate,
            nitrogen_half_saturation,
            alpha,
        )
    )
end

"""
Build expression for a simplified detritus function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function typical_detritus(plankton_array)
    return :(
        net_linear_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            linear_mortality,
            1 - mortality_export_fraction,
        ) +
        net_preferential_predation_assimilation_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            holling_half_saturation,
            maximum_predation_rate,
            assimilation_efficiency_matrix,
            palatability_matrix,
        ) +
        net_quadratic_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            quadratic_mortality,
            1 - mortality_export_fraction,
        ) - idealized_remineralization(D, detritus_remineralization)
    )
end

"""
Build expression for a simplified phytoplankton growth function.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
- `plankton_name`: name of the phytoplankton for which we are returning the expression passed
    as a String (e.g., "P1").
"""
function simplified_phytoplankton_growth(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        idealized_photosynthetic_growth(
            N,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_name],
            nitrogen_half_saturation[$plankton_name],
            alpha[$plankton_name],
        ) - summed_preferential_predation_loss(
            $plankton_name,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            maximum_predation_rate,
            holling_half_saturation,
            palatability_matrix,
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_name])
    )
end

"""
Build expression for simplified zooplankton growth function.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
- `plankton_name`: name of the zooplankton for which we are returning the expression passed
    as a String (e.g., "Z1").
"""
function simplified_zooplankton_growth(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        summed_preferential_predation_gain(
            $plankton_name,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            assimilation_efficiency_matrix,
            maximum_predation_rate,
            holling_half_saturation,
            palatability_matrix,
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_name]) -
        quadratic_loss($(plankton_symbol), quadratic_mortality[$plankton_name])
    )
end

end # module
