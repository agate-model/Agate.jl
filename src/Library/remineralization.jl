module Remineralization

export idealized_remineralization

"""
Idealized remineralization of detritus into nutrients.

# Arguments
- `D`: detritus
- `r`: remineralization rate
"""
function idealized_remineralization(D, r)
    return r * D
end

end # module
