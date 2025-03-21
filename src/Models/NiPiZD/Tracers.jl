module Tracers

using NamedArrays

export detritus_default,
    nutrients_default,
    nutrients_geider_light,
    phytoplankton_default,
    phytoplankton_geider_light,
    zooplankton_default

"""
Build expression for a single nutrient function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function nutrients_default(plankton_array)
    return :(
        sum(
            linear_loss.([$(plankton_array...)], linear_mortality.array) *
            mortality_export_fraction,
        ) +
        sum(
            quadratic_loss.([$(plankton_array...)], quadratic_mortality.array) *
            mortality_export_fraction,
        ) +
        remineralization_idealized(D, detritus_remineralization) - sum(
            photosynthetic_growth_single_nutrient.(
                N,
                [$(plankton_array...)],
                PAR,
                # extract an array of values from the NamedArray
                # otherwise get an error when trying to broadcast
                maximum_growth_rate.array,
                nutrient_half_saturation.array,
                alpha.array,
            ),
        )
    )
end

"""
Build expression for a single nutrient function of time where photosynthetic growth is limited based on the Geider formulation.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function nutrients_geider_light(plankton_array)
    return :(
        sum(
            linear_loss.([$(plankton_array...)], linear_mortality.array) *
            mortality_export_fraction,
        ) +
        sum(
            quadratic_loss.([$(plankton_array...)], quadratic_mortality.array) *
            mortality_export_fraction,
        ) +
        remineralization_idealized(D, detritus_remineralization) - sum(
            photosynthetic_growth_single_nutrient_geider_light.(
                N,
                [$(plankton_array...)],
                PAR,
                # size dependant values
                maximum_growth_rate.array,
                nutrient_half_saturation.array,
                photosynthetic_slope.array,
                chlorophyll_to_carbon_ratio.array,
            ),
        )
    )
end

"""
Build expression for a simplified detritus function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `phyto_array`: names of all phytoplankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2]`
- `zoo_array`: names of all zooplankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2]`
"""
function detritus_default(phyto_array, zoo_array)
    plankton_array = vcat(zoo_array, phyto_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality["P1"])) *
        (1 - mortality_export_fraction) +
        sum(linear_loss.([$(zoo_array...)], linear_mortality["Z1"])) *
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
                assimilation_efficiency_matrix.array,
                # predator size dependant parameters -> apply column-wise
                maximum_predation_rate.array,
                holling_half_saturation["Z1"],
                # predator x prey matrix
                palatability_matrix.array,
            ),
        ) +
        sum(
            quadratic_loss.([$(zoo_array...)], quadratic_mortality["Z1"]) *
            (1 - mortality_export_fraction),
        ) - remineralization_idealized(D, detritus_remineralization)
    )
end

"""
Build expression for a simplified phytoplankton growth function.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
- `plankton_name`: name of the phytoplankton for which we are returning the expression passed
    as a String (e.g., "P1").
"""
function phytoplankton_default(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        photosynthetic_growth_single_nutrient(
            N,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_name],
            nutrient_half_saturation[$plankton_name],
            alpha["P1"],
        ) - sum(
            predation_loss_preferential.(
                # the prey
                $(plankton_symbol),
                # all potential predators
                [$(plankton_array...)],
                # predator size dependant parameters
                maximum_predation_rate.array,
                holling_half_saturation["Z1"],
                # get the prey column -> sum over all predator rows
                palatability_matrix[:, $plankton_name].array,
            ),
        ) - linear_loss($(plankton_symbol), linear_mortality["P1"])
    )
end

"""
Build expression for a simplified phytoplankton growth function which adds geider light limitation.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
- `plankton_name`: name of the phytoplankton for which we are returning the expression passed
    as a String (e.g., "P1").
"""
function phytoplankton_geider_light(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        photosynthetic_growth_single_nutrient_geider_light(
            N,
            $(plankton_symbol),
            PAR,
            maximum_growth_rate[$plankton_name],
            nutrient_half_saturation[$plankton_name],
            photosynthetic_slope["P1"],
            chlorophyll_to_carbon_ratio["P1"],
        ) - sum(
            # exactly the same as in phytoplankton_growth_single_nutrient
            predation_loss_preferential.(
                $(plankton_symbol),
                [$(plankton_array...)],
                maximum_predation_rate.array,
                holling_half_saturation["Z1"],
                palatability_matrix[:, $plankton_name].array,
            ),
        ) - linear_loss($(plankton_symbol), linear_mortality["P1"])
    )
end

"""
Build expression for simplified zooplankton growth function.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
- `plankton_name`: name of the zooplankton for which we are returning the expression passed
    as a String (e.g., "Z1").
"""
function zooplankton_default(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        sum(
            predation_gain_preferential.(
                # all potential prey
                [$(plankton_array...)],
                # the predator
                $(plankton_symbol),
                # get the predator row -> sum over all prey columns
                assimilation_efficiency_matrix[$plankton_name, :],
                # predator size dependant parameter
                maximum_predation_rate[$plankton_name],
                holling_half_saturation["Z1"],
                # get the predator row -> sum over all prey columns
                palatability_matrix[$plankton_name, :],
            ),
        ) - linear_loss($(plankton_symbol), linear_mortality["Z1"]) -
        quadratic_loss($(plankton_symbol), quadratic_mortality["Z1"])
    )
end

end # module
