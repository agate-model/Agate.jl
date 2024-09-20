"""
Modules related to plankton nutrient uptake

"""
module Nutrients

export
    monod_limitation

"
    N / (kₙ + N)

Monod formulation of nutrient limitation, which is based on Michaelis-Menten enzyme kinetics. 

Where: 
R = nutrient (e.g. N, P, Si)
kᵣ = nutrient half saturation constant

Note that sometimes this formulation is also used for Predation.

"
function monod_limitation(R, kᵣ)
    R / (kᵣ + R)
end

end # module