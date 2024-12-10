module Tracers

using NamedArrays

export typical_detritus,
    typical_nutrients, simplified_phytoplankton_growth, simplified_zooplankton_growth

"""
Build expression for a single nutrient....

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
        net_photosynthetic_growth(
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
Build expression for a simplified detritus...

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
        net_predation_assimilation_loss(
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
Build expression for a simplified phytoplankton growth function....

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`

# Notes on args in the expression
- `plankton_name`: The name of the plankton for which the rate of change is estimated
- `N`: Nitrogen
- `P`: NamedArray which includes all plankton
- `linear_mortality`: NamedArray of all plankton linear mortality rates
- `quadratic_mortality`: ....
- `maximum_growth_rate`: NamedArray of all plankton maximum growth rates
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `maximum_predation_rate`: NamedArray of all plankton maximum predation rates
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
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
        ) - summed_predation_loss(
            $plankton_name,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            maximum_predation_rate,
            holling_half_saturation,
            palatability_matrix,
        ) - linear_loss($(plankton_symbol), linear_mortality[$plankton_name]) -
        quadratic_loss($(plankton_symbol), quadratic_mortality[$plankton_name])
    )
end

"""
Build expression for simplified zooplankton growth function...

# Arguments
- `plankton_array`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:P1, :P2, :Z1, :Z2]`

# Notes on args in the expression
- `plankton_name`: The name of the plankton for which the rate of change is estimated
- `P`: NamedArray which includes all plankton
- `linear_mortality`: NamedArray of all plankton linear mortality rates
- `quadratic_mortality`: ....
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `maximum_predation_rate`: NamedArray of all plankton maximum predation rates
- `assimilation efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function simplified_zooplankton_growth(plankton_array, plankton_name)
    plankton_symbol = Symbol(plankton_name)
    return :(
        summed_predation_gain(
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
