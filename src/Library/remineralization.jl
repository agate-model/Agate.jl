"""Remineralization kernels."""

module Remineralization

export linear_remineralization

"""
    linear_remineralization(D, rate)

Idealized remineralization of detritus into dissolved nutrients.

!!! formulation
    r * D

    where:
    - D = detritus concentration
    - r = remineralization rate

# Arguments
- `D`: detritus concentration
- `rate`: remineralization rate
"""
@inline linear_remineralization(D, rate) = rate * D

end # module
