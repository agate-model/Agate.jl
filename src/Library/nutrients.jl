"""
Functions related to plankton nutrient uptake.
"""

module Nutrients

export monod_limitation

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

end # module
