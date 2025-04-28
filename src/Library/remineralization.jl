module Remineralization

export remineralization_idealized

"""
    remineralization_idealized(D, remineralization_rate)
    
Idealized remineralization of detritus into dissolved nutrients.

!!! formulation
    r * D

    where:
    - D = detritus concentration
    - r = remineralization rate

# Arguments
- `D`: detritus concentration
- `remineralization_rate`: remineralization rate
"""
function remineralization_idealized(D, remineralization_rate)
    return remineralization_rate * D
end

end # module
