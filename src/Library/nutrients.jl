"""Functions related to plankton nutrient uptake."""

module Nutrients

export monod_limitation, liebig_minimum

"""
    monod_limitation(nutrient_concentration, nutrient_half_saturation)

Monod formulation of nutrient limitation:

```math
R / (k + R)
```
"""
@inline function monod_limitation(nutrient_concentration, nutrient_half_saturation)
    return nutrient_concentration / (nutrient_half_saturation + nutrient_concentration)
end

"""
    liebig_minimum(a, b, rest...)

Liebig's law of the minimum.

This implementation provides allocation-free overloads for scalar inputs.
"""
@inline liebig_minimum(a, b) = ifelse(a < b, a, b)

@inline function liebig_minimum(a, b, c, rest...)
    return liebig_minimum(liebig_minimum(a, b), c, rest...)
end

@inline function liebig_minimum(values::NTuple{N, T}) where {N, T}
    m = values[1]
    @inbounds for i in 2:N
        m = liebig_minimum(m, values[i])
    end
    return m
end

end # module
