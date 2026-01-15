module Mortality

using Agate.Equations: AExpr, sum_over
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

# -----------------------------------------------------------------------------
# Symbolic (construction-time) equation blocks
# -----------------------------------------------------------------------------

"""\
    linear_loss(PV, plankton_sym::Symbol, idx::Int) -> AExpr

Symbolic linear mortality loss term for a single plankton size-class.

Missing/`nothing` behavior is controlled by the model's parameter registry.
"""
function linear_loss(PV, plankton_sym::Symbol, idx::Int)
    return PV.linear_mortality[idx] * plankton_sym
end

"""\
    quadratic_loss(PV, plankton_sym::Symbol, idx::Int) -> AExpr

Symbolic quadratic mortality loss term for a single plankton size-class.

Missing/`nothing` behavior is controlled by the model's parameter registry.
"""
function quadratic_loss(PV, plankton_sym::Symbol, idx::Int)
    P = plankton_sym
    # Avoid evaluating `P * P` at construction time (which would try to multiply Symbols).
    # Left-associative `AExpr * Symbol * Symbol` builds a symbolic Expr safely.
    return PV.quadratic_mortality[idx] * P * P
end


"""
    linear_loss_sum(PV, plankton_syms)

Build an allocation-free `Expr` summing linear mortality over all plankton symbols.

The expression expects the runtime container `linear_mortality` to be in scope (a vector of length `n_total`).
"""
function linear_loss_sum(PV, plankton_syms::AbstractVector{Symbol})
    return sum_over(plankton_syms) do sym, i
        PV.linear_mortality[i] * sym
    end
end

"""
    quadratic_loss_sum(PV, plankton_syms)

Build an allocation-free `Expr` summing quadratic mortality over all plankton symbols.

The expression expects the runtime container `quadratic_mortality` to be in scope (a vector of length `n_total`).
"""
function quadratic_loss_sum(PV, plankton_syms::AbstractVector{Symbol})
    return sum_over(plankton_syms) do sym, i
        # Avoid evaluating `sym * sym` at construction time; build it symbolically.
        PV.quadratic_mortality[i] * sym * sym
    end
end

end # module
