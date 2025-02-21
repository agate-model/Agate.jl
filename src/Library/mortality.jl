module Mortality

export linear_loss, quadratic_loss, net_linear_loss, net_quadratic_loss

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

"""
Net loss of all plankton due to linear mortality.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `linear_mortality`: NamedArray of plankton linear mortality rates
- `fraction`:
"""
function net_linear_loss(P, linear_mortality, fraction)
    return sum(linear_mortality * P) * fraction
end

"""
Net loss of all plankton due to quadratic mortality.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `quadratic_mortality`: plankton quadratic mortality rate
- `fraction`: `
"""
function net_quadratic_loss(P, quadratic_mortality, fraction)
    return sum(quadratic_mortality * (P .^ 2)) * fraction
end

end # module
