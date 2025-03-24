module Tracers

export detritus_default,
    nutrients_default,
    nutrients_geider_light,
    phytoplankton_default,
    phytoplankton_geider_light,
    zooplankton_default

"""
    nutrients_default(plankton_array)

Build expression for a single nutrient function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function nutrients_default(plankton_array)
    return :(
        sum(
            linear_loss.([$(plankton_array...)], linear_mortality) *
            mortality_export_fraction,
        ) +
        sum(
            quadratic_loss.([$(plankton_array...)], quadratic_mortality) *
            mortality_export_fraction,
        ) +
        remineralization_idealized(D, detritus_remineralization) - sum(
            photosynthetic_growth_single_nutrient.(
                N,
                [$(plankton_array...)],
                PAR,
                maximum_growth_rate,
                nutrient_half_saturation,
                alpha,
            ),
        )
    )
end

"""
    nutrients_geider_light(plankton_array)

Build expression for a single nutrient function of time where photosynthetic growth is limited
based on the Geider formulation.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function nutrients_geider_light(plankton_array)
    return :(
        sum(
            linear_loss.([$(plankton_array...)], linear_mortality) *
            mortality_export_fraction,
        ) +
        sum(
            quadratic_loss.([$(plankton_array...)], quadratic_mortality) *
            mortality_export_fraction,
        ) +
        remineralization_idealized(D, detritus_remineralization) - sum(
            photosynthetic_growth_single_nutrient_geider_light.(
                N,
                [$(plankton_array...)],
                PAR,
                # size dependant values
                maximum_growth_rate,
                nutrient_half_saturation,
                photosynthetic_slope,
                chlorophyll_to_carbon_ratio,
            ),
        )
    )
end

"""
    detritus_default(plankton_array)

Build expression for a simplified detritus function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function detritus_default(plankton_array)
    return :(
        sum(linear_loss.([$(plankton_array...)], linear_mortality)) *
        (1 - mortality_export_fraction) +
        sum(
            # the function includes predator x prey matrix inputs so have to make sure that
            # arrays are broadcast row/column wise as appropriate
            predation_assimilation_loss_preferential.(
                # prey is row-wise -> use array' to declare a row vector
                [$(plankton_array...)]',
                # predators are column-wise so leave as is
                [$(plankton_array...)],
                # predator x prey matrix
                assimilation_efficiency_matrix,
                # predator size dependant parameters -> apply column-wise
                maximum_predation_rate,
                holling_half_saturation,
                # predator x prey matrix
                palatability_matrix,
            ),
        ) +
        sum(
            quadratic_loss.([$(plankton_array...)], quadratic_mortality) *
            (1 - mortality_export_fraction),
        ) - remineralization_idealized(D, detritus_remineralization)
    )
end

"""
    phytoplankton_default(plankton_array, plankton_name, plankton_idx)

Build expression for a simplified phytoplankton growth function.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
- `plankton_name`: name of the phytoplankton for which we are returning the expression, passed
    as a String (e.g., "P1").
- `plankton_idx`: the index at which this plankton's values are stored in all parameter Arrays
"""
function phytoplankton_default(plankton_array, plankton_name, plankton_idx)
    plankton_symbol = Symbol(plankton_name)
    return :(
        photosynthetic_growth_single_nutrient(
            N,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_idx],
            nutrient_half_saturation[$plankton_idx],
            alpha[$plankton_idx],
        ) - sum(
            predation_loss_preferential.(
                # the prey
                $(plankton_symbol),
                # all potential predators
                [$(plankton_array...)],
                # predator size dependant parameters
                maximum_predation_rate,
                holling_half_saturation,
                # get the prey column -> sum over all predator rows
                palatability_matrix[:, $plankton_idx],
            ),
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_idx])
    )
end

"""
    phytoplankton_geider_light(plankton_array, plankton_name, plankton_idx)

Build expression for a simplified phytoplankton growth function which adds geider light limitation.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
- `plankton_name`: name of the phytoplankton for which we are returning the expression, passed
    as a String (e.g., "P1").
- `plankton_idx`: the index at which this plankton's values are stored in all parameter Arrays
"""
function phytoplankton_geider_light(plankton_array, plankton_name, plankton_idx)
    plankton_symbol = Symbol(plankton_name)
    return :(
        photosynthetic_growth_single_nutrient_geider_light(
            N,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_idx],
            nutrient_half_saturation[$plankton_idx],
            photosynthetic_slope[$plankton_idx],
            chlorophyll_to_carbon_ratio[$plankton_idx],
        ) - sum(
            # exactly the same as in phytoplankton_growth_single_nutrient
            predation_loss_preferential.(
                $(plankton_symbol),
                [$(plankton_array...)],
                maximum_predation_rate,
                holling_half_saturation,
                palatability_matrix[:, $plankton_idx],
            ),
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_idx])
    )
end

"""
    zooplankton_default(plankton_array, plankton_name, plankton_idx)

Build expression for simplified zooplankton growth function.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
- `plankton_name`: name of the zooplankton for which we are returning the expression passed
    as a String (e.g., "Z1").
- `plankton_idx`: the index at which this plankton's values are stored in all parameter Arrays
"""
function zooplankton_default(plankton_array, plankton_name, plankton_idx)
    plankton_symbol = Symbol(plankton_name)
    return :(
        sum(
            predation_gain_preferential.(
                # all potential prey
                [$(plankton_array...)],
                # the predator
                $(plankton_symbol),
                # get the predator row -> sum over all prey columns
                assimilation_efficiency_matrix[$plankton_idx, :],
                # predator size dependant parameter
                maximum_predation_rate[$plankton_idx],
                holling_half_saturation[$plankton_idx],
                # get the predator row -> sum over all prey columns
                palatability_matrix[$plankton_idx, :],
            ),
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_idx]) -
        quadratic_loss($(plankton_symbol), quadratic_mortality[$plankton_idx])
    )
end

end # module
