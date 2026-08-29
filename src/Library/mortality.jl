"""Mortality and loss kernels."""

module Mortality

export linear_loss

"""
    linear_loss(P, rate)

Linear mortality (loss) rate.

!!! formulation
    ``l`` * ``P``

    where:
    - ``P`` = plankton concentration
    - ``l`` = mortality (loss) rate

# Arguments
- `P`: plankton concentration
- `rate`: mortality (loss) rate
"""
@inline linear_loss(P, rate) = rate * P

end # module
