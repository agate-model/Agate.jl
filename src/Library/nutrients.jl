"""
Functions related to plankton nutrient uptake.
"""

module Nutrients

export monod_limitation, liebig_minimum

"
    R / (kᵣ + R)

Monod formulation of nutrient limitation, which is based on Michaelis-Menten enzyme kinetics.

# Arguments
- `R`: nutrient (e.g. N, P, Si)
- `kᵣ`: nutrient half saturation constant

Note that sometimes this formulation is also used for Predation ('Holling type 2').
"
function monod_limitation(R, kᵣ)
    return R / (kᵣ + R)
end

"""
    liebig_minimum(nutrient_limitations)

Liebig's law of the minimum, which states that growth is limited by the scarcest (most limiting) resource.

# Arguments
- `nutrient_limitations`: an array of nutrient limitation values

Returns the minimum value among the given nutrient limitations.
"""
function liebig_minimum(nutrient_limitations)
    return minimum(nutrient_limitations)
end

end # module
