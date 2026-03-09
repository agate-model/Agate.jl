"""Remineralization functors."""

module Remineralization

export linear_remineralization

"""
    LinearRemineralization(rate)

Linear remineralization functor.

!!! formulation
    r * D

    where:
    - D = detritus (or organic matter) concentration
    - r = remineralization rate
"""
struct LinearRemineralization{T}
    rate::T
end

@inline (r::LinearRemineralization)(D) = r.rate * D

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
@inline linear_remineralization(D, rate) = LinearRemineralization(rate)(D)

end # module
