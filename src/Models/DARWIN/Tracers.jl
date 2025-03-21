module Tracers

using NamedArrays

export detritus_typical,
    DIC_geider_light,
    DIN_geider_light,
    PO4_geider_light,
    POC_default,
    DOC_default,
    PON_default,
    DON_default,
    POP_default,
    DOP_default,
    phytoplankton_growth_two_nutrients_geider_light,
    zooplankton_default

"""
    DIC_geider_light(phyto_array)

Build expression representing the evolution of DIC over time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `phyto_array`: names of all the phytoplankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2]`
"""
function DIC_geider_light(phyto_array)
    return :(
        remineralization_idealized(DOC, DOC_remineralization) +
        remineralization_idealized(POC, POC_remineralization) - sum(
            photosynthetic_growth_two_nutrients_geider_light.(
                DIN,
                PO4,
                [$(phyto_array...)],
                PAR,
                maximum_growth_rate.array,
                half_saturation_DIN.array,
                half_saturation_PO4.array,
                photosynthetic_slope,
                chlorophyll_to_carbon_ratio,
            ),
        )
    )
end

"""
    DIN_geider_light(phyto_array)

Build expression representing the evolution of DIN over time assuming fixed stoichiometry.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `phyto_array`: names of all the phytoplankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2]`
"""
function DIN_geider_light(phyto_array)
    return :(
        remineralization_idealized(DOC, DOC_remineralization) +
        remineralization_idealized(POC, POC_remineralization) - sum(
            photosynthetic_growth_two_nutrients_geider_light.(
                DIN,
                PO4,
                [$(phyto_array...)],
                PAR,
                maximum_growth_rate.array,
                half_saturation_DIN.array,
                half_saturation_PO4.array,
                photosynthetic_slope,
                chlorophyll_to_carbon_ratio,
            ) * nitrogen_to_carbon,
        )
    )
end

"""
    PO4_geider_light(phyto_array)

Build expression representing the evolution of DIN over time assuming fixed stoichiometry.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `phyto_array`: names of all the phytoplankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2]`
"""
function PO4_geider_light(phyto_array)
    return :(
        remineralization_idealized(DOC, DOC_remineralization) +
        remineralization_idealized(POC, POC_remineralization) - sum(
            photosynthetic_growth_two_nutrients_geider_light.(
                DIN,
                PO4,
                [$(phyto_array...)],
                PAR,
                maximum_growth_rate.array,
                half_saturation_DIN.array,
                half_saturation_PO4.array,
                photosynthetic_slope,
                chlorophyll_to_carbon_ratio,
            ) * nitrogen_to_carbon,
        )
    )
end

"""
    DOC_default(phyto_array, zoo_array)

Build expression for a simplified DOC function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `phyto_array`: names of all phytoplankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2]`
- `zoo_array`: names of all zooplankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2]`
"""
function DOC_default(phyto_array, zoo_array)
    plankton_array = vcat(zoo_array, phyto_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality_p)) *
        (1 - DOM_POM_fractionation) +
        sum(linear_loss.([$(zoo_array...)], linear_mortality_z)) *
        (1 - DOM_POM_fractionation) +
        sum(
            # essentially same as the detritus_typical function
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
                holling_half_saturation,
                # predator x prey matrix
                palatability_matrix.array,
            ),
        ) * (1 - DOM_POM_fractionation) +
        sum(
            quadratic_loss.([$(zoo_array...)], quadratic_mortality) *
            (1 - DOM_POM_fractionation),
        ) - remineralization_idealized(DOC, DOC_remineralization)
    )
end

"""
    DON_default(phyto_array, zoo_array)

Build expression for a simplified DON function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `phyto_array`: names of all phytoplankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2]`
- `zoo_array`: names of all zooplankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2]`
"""
function DON_default(phyto_array, zoo_array)
    plankton_array = vcat(zoo_array, phyto_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality_p)) *
        (1 - DOM_POM_fractionation) *
        nitrogen_to_carbon +
        sum(linear_loss.([$(zoo_array...)], linear_mortality_z)) *
        (1 - DOM_POM_fractionation) *
        nitrogen_to_carbon +
        sum(
            # essentially same as the detritus_typical function
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
                holling_half_saturation,
                # predator x prey matrix
                palatability_matrix.array,
            ) * nitrogen_to_carbon,
        ) * (1 - DOM_POM_fractionation) +
        sum(
            quadratic_loss.([$(zoo_array...)], quadratic_mortality) *
            (1 - DOM_POM_fractionation) *
            nitrogen_to_carbon,
        ) - remineralization_idealized(DON, DON_remineralization)
    )
end

"""
    DOP_default(phyto_array, zoo_array)

Build expression for a simplified DOP function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `phyto_array`: names of all phytoplankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2]`
- `zoo_array`: names of all zooplankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2]`
"""
function DOP_default(phyto_array, zoo_array)
    plankton_array = vcat(zoo_array, phyto_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality_p)) *
        (1 - DOM_POM_fractionation) *
        phosphorus_to_carbon +
        sum(linear_loss.([$(zoo_array...)], linear_mortality_z)) *
        (1 - DOM_POM_fractionation) *
        phosphorus_to_carbon +
        sum(
            # essentially same as the detritus_typical function
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
                holling_half_saturation,
                # predator x prey matrix
                palatability_matrix.array,
            ) * phosphorus_to_carbon,
        ) * (1 - DOM_POM_fractionation) +
        sum(
            quadratic_loss.([$(zoo_array...)], quadratic_mortality) *
            (1 - DOM_POM_fractionation) *
            phosphorus_to_carbon,
        ) - remineralization_idealized(DOP, DOP_remineralization)
    )
end

"""
    POC_default(phyto_array, zoo_array)

Build expression for a simplified POC function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `phyto_array`: names of all phytoplankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2]`
- `zoo_array`: names of all zooplankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2]`
"""
function POC_default(phyto_array, zoo_array)
    plankton_array = vcat(zoo_array, phyto_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality_p)) * DOM_POM_fractionation +
        sum(linear_loss.([$(zoo_array...)], linear_mortality_z)) * DOM_POM_fractionation +
        sum(
            # essentially same as the detritus_typical function
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
                holling_half_saturation,
                # predator x prey matrix
                palatability_matrix.array,
            ),
        ) * DOM_POM_fractionation +
        sum(quadratic_loss.([$(zoo_array...)], quadratic_mortality) * DOM_POM_fractionation) -
        remineralization_idealized(POC, POC_remineralization)
    )
end

"""
    PON_default(phyto_array, zoo_array)

Build expression for a simplified PON function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`
"""
function PON_default(phyto_array, zoo_array)
    plankton_array = vcat(zoo_array, phyto_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality_p)) *
        DOM_POM_fractionation *
        nitrogen_to_carbon +
        sum(linear_loss.([$(zoo_array...)], linear_mortality_z)) *
        DOM_POM_fractionation *
        nitrogen_to_carbon +
        sum(
            # essentially same as the detritus_typical function
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
                holling_half_saturation,
                # predator x prey matrix
                palatability_matrix.array,
            ),
        ) *
        DOM_POM_fractionation *
        nitrogen_to_carbon +
        sum(
            quadratic_loss.([$(zoo_array...)], quadratic_mortality) *
            DOM_POM_fractionation *
            nitrogen_to_carbon,
        ) - remineralization_idealized(PON, PON_remineralization)
    )
end

"""
    POP_default(phyto_array, zoo_array)

Build expression for a simplified POP function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for vectorization to work (and arranged in the same plankton order).

# Arguments
- `phyto_array`: names of all phytoplankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2]`
- `zoo_array`: names of all zooplankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2]`
"""
function POP_default(phyto_array, zoo_array)
    plankton_array = vcat(zoo_array, phyto_array)
    return :(
        sum(linear_loss.([$(phyto_array...)], linear_mortality_p)) *
        DOM_POM_fractionation *
        phosphorus_to_carbon +
        sum(linear_loss.([$(zoo_array...)], linear_mortality_z)) *
        DOM_POM_fractionation *
        phosphorus_to_carbon +
        sum(
            # essentially same as the detritus_typical function
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
                holling_half_saturation,
                # predator x prey matrix
                palatability_matrix.array,
            ),
        ) *
        DOM_POM_fractionation *
        phosphorus_to_carbon +
        sum(
            quadratic_loss.([$(zoo_array...)], quadratic_mortality) *
            DOM_POM_fractionation *
            phosphorus_to_carbon,
        ) - remineralization_idealized(POP, POP_remineralization)
    )
end

"""
    phytoplankton_growth_two_nutrients_geider_light(plankton_array, plankton_name)

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
        ) - sum(
            predation_loss_preferential.(
                # the prey
                $(plankton_symbol),
                # all potential predators
                [$(plankton_array...)],
                # predator size dependant parameters
                maximum_predation_rate.array,
                holling_half_saturation.array,
                # get the prey column -> sum over all predator rows
                palatability_matrix[:, $plankton_name].array,
            ),
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_name])
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
                holling_half_saturation[$plankton_name],
                # get the predator row -> sum over all prey columns
                palatability_matrix[$plankton_name, :],
            ),
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_name]) -
        quadratic_loss($(plankton_symbol), quadratic_mortality[$plankton_name])
    )
end

end # module
