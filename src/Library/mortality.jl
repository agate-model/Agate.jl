module Mortality

using ..ExprUtils: sum_expr

export linear_loss, quadratic_loss, linear_loss_sum, quadratic_loss_sum

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


"""
    linear_loss_sum(plankton_syms)

Build an allocation-free `Expr` summing linear mortality over all plankton symbols.

The expression expects the runtime container `linear_mortality` to be in scope (a vector of length `n_total`).
"""
function linear_loss_sum(plankton_syms::AbstractVector{Symbol})
    terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(terms, :(linear_loss($sym, linear_mortality[$i])))
    end
    return sum_expr(terms)
end

"""
    quadratic_loss_sum(plankton_syms)

Build an allocation-free `Expr` summing quadratic mortality over all plankton symbols.

The expression expects the runtime container `quadratic_mortality` to be in scope (a vector of length `n_total`).
"""
function quadratic_loss_sum(plankton_syms::AbstractVector{Symbol})
    terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(terms, :(quadratic_loss($sym, quadratic_mortality[$i])))
    end
    return sum_expr(terms)
end

end # module
