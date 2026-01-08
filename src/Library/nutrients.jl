"""Functions related to plankton nutrient uptake."""

module Nutrients

export monod_limitation, liebig_minimum

"""
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

!!! tip
    Sometimes this formulation is also used for predation (≈'Holling type 2').
"""
@inline function monod_limitation(nutrient_concentration, nutrient_half_saturation)
    return nutrient_concentration / (nutrient_half_saturation + nutrient_concentration)
end

"""
    liebig_minimum(a, b, rest...)

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
@inline liebig_minimum(a, b) = ifelse(a < b, a, b)

@inline function liebig_minimum(a, b, c, rest...)
    return liebig_minimum(liebig_minimum(a, b), c, rest...)
end

@inline function liebig_minimum(values::NTuple{N,T}) where {N,T}
    m = values[1]
    @inbounds for i in 2:N
        m = liebig_minimum(m, values[i])
    end
    return m
end

end # module
