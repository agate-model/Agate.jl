module Mortality

export linear_loss, quadratic_loss, net_linear_loss, net_quadratic_loss, net_linear_loss_quota, net_quadratic_loss_quota

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
- `linear_mortality`: NamedArray of all plankton linear mortality rates
"""
function net_linear_loss(P, linear_mortality, DOM_POM_fractionation)
    # sum over all plankton in `P`
    print(DOM_POM_fractionation)
    return sum([linear_loss(P[name], linear_mortality[name])*DOM_POM_fractionation for name in names(linear_mortality, 1)])
end

"""
Net loss of all plankton due to linear mortality with a quota term.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `linear_mortality`: NamedArray of all plankton linear mortality rates
- `quota`: NamedArray which includes all plankton quota values
"""
function net_linear_loss_quota(P, linear_mortality, DOM_POM_fractionation, quota)
    # sum over all plankton in `P`
    return sum([linear_loss(P[name], linear_mortality[name])*quota[name]*DOM_POM_fractionation for name in names(P, 1)]) 
end

"""
Net loss of all plankton due to quadratic mortality.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `quadratic_mortality`: NamedArray of all plankton quadratic mortality rates
"""
function net_quadratic_loss(P, quadratic_mortality, DOM_POM_fractionation)
    # sum over plankton that have a `quadratic_mortality`
    return sum(
        [
            quadratic_loss(P[name], quadratic_mortality[name]) * DOM_POM_fractionation for
            name in names(quadratic_mortality, 1)
        ] 
    )
end

"""
Net loss of all plankton due to quadratic mortality with a quota term.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `quadratic_mortality`: NamedArray of all plankton quadratic mortality rates
- `quota`: NamedArray which includes all plankton quota values
"""
function net_quadratic_loss_quota(P, quadratic_mortality, DOM_POM_fractionation, quota)
    # sum over plankton that have a `quadratic_mortality`
    return sum(
        [
            quadratic_loss(P[name], quadratic_mortality[name])*quota[name]*DOM_POM_fractionation for  
            name in names(quadratic_mortality, 1)
        ] 
    )
end

end # module
