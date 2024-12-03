module Mortality

export linear_loss, quadratic_loss

"""
Linear mortality rate. 
In this formulation mortality is constant, and can be interpreted as 
a "closure term" for low density predation and and other death terms.

# Arguments
- `P`: plankton concentration
- `l`: mortality rate
"""
linear_loss(P, l) = l * P

"""
Quadratic mortality coefficient.
In this formulation mortality increases exponentially with plankton biomass
and is often interpreted to represent viral processes and non-represented density-dependent predation effects.

# Arguments
- `P`: plankton concentration
- `l`: mortality rate
"""
quadratic_loss(P, l) = l * P^2

end # module
