module Mortality

export linear_loss, quadratic_loss

"""
    linear_loss(plankton_concentration, mortality)

Linear mortality rate.

!!! formulation
    ``l`` * ``P``

    where:
    - ``P`` = plankton concentration
    - ``l`` = mortality rate

In this formulation mortality is constant, and can be interpreted as
a "closure term" for low density predation and and other death terms.

# Arguments
- `plankton_concentration`: plankton concentration
- `mortality`: mortality rate
"""
linear_loss(plankton_concentration, mortality) = mortality * plankton_concentration

"""
    quadratic_loss(plankton_concentration, mortality)

Quadratic mortality rate.

!!! formulation

    ``l`` * ``P``^2

    where:
    - ``P`` = plankton concentration
    - ``l`` = mortality rate

In this formulation mortality increases exponentially with plankton biomass
and is often interpreted to represent viral processes and non-represented density-dependent predation effects.

# Arguments
- `plankton_concentration`: plankton concentration
- `mortality`: mortality rate
"""
quadratic_loss(plankton_concentration, mortality) = mortality * plankton_concentration^2

end # module
