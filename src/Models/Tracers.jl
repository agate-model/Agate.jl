module Tracers

using NamedArrays

export detritus_typical,
    nutrients_typical,
    nutrients_geider_light,
    phytoplankton_growth_single_nutrient,
    phytoplankton_growth_single_nutrient_geider_light,
    zooplankton_growth_simplified

"""
Build expression for a single nutrient function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function nutrients_typical(phyto_array, zoo_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality["P"]) * mortality_export_fraction) +
        sum(linear_loss.([$(zoo_array...)], linear_mortality["Z"]) * mortality_export_fraction)  +
        sum(quadratic_loss.([$(zoo_array...)], quadratic_mortality["Z"]) * mortality_export_fraction) +
        remineralization_idealized(D, detritus_remineralization) -
        sum(photosynthetic_growth_single_nutrient.(
            N,
            [$(phyto_array...)],
            PAR,
            maximum_growth_rate.array,
            nutrient_half_saturation.array,
            alpha["P"],
        ))
    )
end

"""
Build expression for a single nutrient function of time where photosynthetic growth is limited based on the Geider formulation.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either a NamedArray or a Float.

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function nutrients_geider_light(phyto_array, zoo_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality["P"]) * mortality_export_fraction) +
        sum(linear_loss.([$(zoo_array...)], linear_mortality["Z"]) * mortality_export_fraction)  +
        sum(quadratic_loss.([$(zoo_array...)], quadratic_mortality["Z"]) * mortality_export_fraction) +
        remineralization_idealized(D, detritus_remineralization) -
        sum(photosynthetic_growth_single_nutrient_geider_light.(
            N,
            [$(phyto_array...)],
            PAR,
            maximum_growth_rate.array,
            nutrient_half_saturation.array,
            photosynthetic_slope["P"],
            chlorophyll_to_carbon_ratio["P"],
        ))
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
function detritus_typical(phyto_array, zoo_array)
    plankton_array = vcat(phyto_array, zoo_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality["P"])) * (1 - mortality_export_fraction) +
        sum(linear_loss.([$(zoo_array...)], linear_mortality["Z"])) * (1 - mortality_export_fraction) +
        net_predation_assimilation_loss_preferential(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            holling_half_saturation,
            maximum_predation_rate,
            assimilation_efficiency_matrix,
            palatability_matrix,
        ) +
        net_quadratic_loss(
            NamedArray([$(zoo_array...)], $(String.(zoo_array))),
            quadratic_mortality["Z"],
            1 - mortality_export_fraction,
        ) - remineralization_idealized(D, detritus_remineralization)
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
function phytoplankton_growth_single_nutrient(phyto_array, zoo_array, plankton_name)
    plankton_array = vcat(phyto_array, zoo_array)
    plankton_symbol = Symbol(plankton_name)
    # remove any digits to get just the type identifier
    plankton_type = replace(plankton_name, r"\d+" => "")
    return :(
        photosynthetic_growth_single_nutrient(
            N,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_name],
            nutrient_half_saturation[$plankton_name],
            alpha[$plankton_type],
        ) - summed_predation_loss_preferential(
            $plankton_name,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            maximum_predation_rate,
            holling_half_saturation,
            palatability_matrix,
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_type])
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
function phytoplankton_growth_single_nutrient_geider_light(phyto_array, zoo_array, plankton_name)
    plankton_array = vcat(phyto_array, zoo_array)
    plankton_symbol = Symbol(plankton_name)
    # remove any digits to get just the type identifier
    plankton_type = replace(plankton_name, r"\d+" => "")
    return :(
        photosynthetic_growth_single_nutrient_geider_light(
            N,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_name],
            nutrient_half_saturation[$plankton_name],
            photosynthetic_slope[$plankton_type],
            chlorophyll_to_carbon_ratio[$plankton_type],
        ) - summed_predation_loss_preferential(
            $plankton_name,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            maximum_predation_rate,
            holling_half_saturation,
            palatability_matrix,
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_type])
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
    # remove any digits to get just the type identifier
    plankton_type = replace(plankton_name, r"\d+" => "")

    return :(
        summed_predation_gain_preferential(
            $plankton_name,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            assimilation_efficiency_matrix,
            maximum_predation_rate,
            holling_half_saturation,
            palatability_matrix,
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_type]) -
        quadratic_loss($(plankton_symbol), quadratic_mortality[$plankton_type])
    )
end

end # module
