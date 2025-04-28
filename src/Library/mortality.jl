module Mortality

export linear_loss, quadratic_loss

"""
    linear_loss(P, mortality)

Linear mortality rate.

!!! formulation
    ``l`` * ``P``

    where:
    - ``P`` = plankton concentration
    - ``l`` = mortality rate

In this formulation mortality is linear, and can be interpreted as
a "closure term" for low density predation and and other death terms.

# Arguments
- `P`: plankton concentration
- `mortality`: mortality rate
"""
linear_loss(P, mortality) = mortality * P

"""
    quadratic_loss(P, mortality)

Quadratic mortality rate.

!!! formulation

    ``l`` * ``P``²

    where:
    - ``P`` = plankton concentration
    - ``l`` = mortality rate

In this formulation mortality increases exponentially with plankton biomass
and is often interpreted to represent viral processes and non-represented density-dependent predation effects.

# Arguments
- `P`: plankton concentration
- `mortality`: mortality rate
"""
quadratic_loss(P, mortality) = mortality * P^2

end # module
