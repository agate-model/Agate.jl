module Mortality

export linear_loss, quadratic_loss

"""
# Arguments
- `P`: plankton concentration
- `l`: mortality rate
"""
linear_loss(P, l) = l * P

"""
# Arguments
- `P`: plankton concentration
- `l`: mortality rate
"""
quadratic_loss(P, l) = l * P^2

end # module
