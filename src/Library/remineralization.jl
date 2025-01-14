module Remineralization

export remineralization_idealized

"""
Idealized remineralization of detritus into dissolved nutrients.

# Arguments
- `D`: detritus
- `r`: remineralization rate
"""
function remineralization_idealized(D, r)
    return r * D
end

end # module
