module Tracers

using NamedArrays

export detritus_typical,
    nutrients_typical,
    nutrients_geider_light,
    DIC_geider_light,
    DIN_geider_light_fixed_ratios,
    PO4_geider_light_fixed_ratios,
    POC_typical,
    DOC_typical,
    phytoplankton_growth_single_nutrient,
    phytoplankton_growth_single_nutrient_geider_light,
    phytoplankton_growth_two_nutrients_geider_light,
    zooplankton_growth_simplified

"""
Build expression for a single nutrient function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function nutrients_typical(plankton_array)
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
        remineralization_idealized(D, detritus_remineralization) -
        net_photosynthetic_growth_single_nutrient(
            N,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            PAR,
            maximum_growth_rate,
            nutrient_half_saturation,
            alpha,
        )
    )
end

"""


The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function DIC_geider_light(plankton_array)
    return :(
        remineralization_idealized(DOC, DOC_remineralization) -
        net_photosynthetic_growth_two_nutrient_geider_light(
            DIN,
            PO4,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            PAR,
            maximum_growth_rate,
            half_saturation_DIN,
            half_saturation_PO4,
            photosynthetic_slope,
            chlorophyll_to_carbon_ratio,
        )
    )
end

"""

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function DIN_geider_light_fixed_ratios(plankton_array)
    return :(
        (
            remineralization_idealized(DOC, DOC_remineralization) -
            net_photosynthetic_growth_two_nutrient_geider_light(
                DIN,
                PO4,
                NamedArray([$(plankton_array...)], $(String.(plankton_array))),
                PAR,
                maximum_growth_rate,
                half_saturation_DIN,
                half_saturation_PO4,
                photosynthetic_slope,
                chlorophyll_to_carbon_ratio,
            )
        ) * nitrogen_to_carbon
    )
end

"""

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function PO4_geider_light_fixed_ratios(plankton_array)
    return :((
        remineralization_idealized(DOC, DOC_remineralization) -
        net_photosynthetic_growth_two_nutrient_geider_light(
            DIN,
            PO4,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            PAR,
            maximum_growth_rate,
            half_saturation_DIN,
            half_saturation_PO4,
            photosynthetic_slope,
            chlorophyll_to_carbon_ratio,
        ) * phosphorus_to_carbon
    ))
end

"""
Build expression for a single nutrient function of time where photosynthetic growth is limited based on the Geider formulation.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function nutrients_geider_light(plankton_array)
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
        remineralization_idealized(D, detritus_remineralization) -
        net_photosynthetic_growth_single_nutrient_geider_light(
            N,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            PAR,
            maximum_growth_rate,
            nutrient_half_saturation,
            photosynthetic_slope,
            chlorophyll_to_carbon_ratio,
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
function detritus_typical(plankton_array)
    return :(
        net_linear_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            linear_mortality,
            1 - mortality_export_fraction,
        ) +
        net_predation_assimilation_loss_preferential(
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
        ) - remineralization_idealized(D, detritus_remineralization)
    )
end

"""
Build expression for a simplified DOC function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function DOC_typical(plankton_array)
    return :(
        net_linear_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            linear_mortality,
            1 - mortality_export_fraction,
        ) +
        net_predation_assimilation_loss_preferential(
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
        ) - remineralization_idealized(DOC, DOC_remineralization)
    )
end

"""
Build expression for a simplified POC function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function POC_typical(plankton_array)
    return :(
        net_linear_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            linear_mortality,
            1 - mortality_export_fraction,
        ) +
        net_predation_assimilation_loss_preferential(
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
        ) - remineralization_idealized(POC, POC_remineralization)
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
function phytoplankton_growth_single_nutrient(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        photosynthetic_growth_single_nutrient(
            N,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_name],
            nutrient_half_saturation[$plankton_name],
            alpha[$plankton_name],
        ) - summed_predation_loss_preferential(
            $plankton_name,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            maximum_predation_rate,
            holling_half_saturation,
            palatability_matrix,
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_name])
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
function phytoplankton_growth_two_nutrients_geider_light(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        photosynthetic_growth_two_nutrients_geider_light(
            DIN,
            PO4,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_name],
            half_saturation_DIN[$plankton_name],
            half_saturation_PO4[$plankton_name],
            photosynthetic_slope[$plankton_name], 
            chlorophyll_to_carbon_ratio[$plankton_name],
        ) - summed_predation_loss_preferential(
            $plankton_name,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            maximum_predation_rate,
            holling_half_saturation,
            palatability_matrix,
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_name])
    )
end

"""
Build expression for a simplified phytoplankton growth function which adds geider light limitation.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
- `plankton_name`: name of the phytoplankton for which we are returning the expression passed
    as a String (e.g., "P1").
"""
function phytoplankton_growth_single_nutrient_geider_light(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        photosynthetic_growth_single_nutrient_geider_light(
            N,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_name],
            nutrient_half_saturation[$plankton_name],
            photosynthetic_slope[$plankton_name],
            chlorophyll_to_carbon_ratio[$plankton_name],
        ) - summed_predation_loss_preferential(
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
function zooplankton_growth_simplified(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        summed_predation_gain_preferential(
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
