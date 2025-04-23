"""
Functions related to plankton nutrient uptake.
"""

module Nutrients

export monod_limitation, liebig_minimum

"
    monod_limitation(nutrient_concentration, nutrient_half_saturation)

Monod formulation of nutrient limitation, which is based on Michaelis-Menten enzyme kinetics.

!!! formulation
    ``R`` / (``kᵣ`` + ``R``)

    where:
    - ``R`` = nutrient concentration (e.g. N, P, Si)
    - ``kᵣ`` = nutrient half saturation constant    

# Arguments
- `nutrient_concentration`: nutrient (e.g. N, P, Si)
- `nutrient_half_saturation`: nutrient half saturation constant

Note that sometimes this formulation is also used for Predation ('Holling type 2').
"
function monod_limitation(nutrient_concentration, nutrient_half_saturation)
    return nutrient_concentration / (nutrient_half_saturation + nutrient_concentration)
end

"""
    liebig_minimum(nutrient_limitations)

Liebig's law of the minimum, which states that growth is limited by the scarcest (most limiting) resource.

!!! formulation
    minimum(nutrient_limitations)

    where:
    - nutrient_limitations = an array of nutrient limitation values
        (e.g. [N, P, Si])

# Arguments
- `nutrient_limitations`: an array of nutrient limitation values

Returns the minimum value among the given nutrient limitations.
"""
function liebig_minimum(nutrient_limitations)
    return minimum(nutrient_limitations)
end

end # module
