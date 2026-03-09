"""Mortality and loss functors."""

module Mortality

export linear_loss, quadratic_loss

"""
    LinearLoss(rate)

Linear mortality / loss functor.

!!! formulation
    ``l`` * ``P``

    where:
    - ``P`` = plankton concentration
    - ``l`` = mortality (loss) rate

In this formulation mortality is linear, and can be interpreted as a closure term for
low-density predation and other linear loss processes.
"""
struct LinearLoss{T}
    rate::T
end

@inline (l::LinearLoss)(P) = l.rate * P

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
@inline linear_loss(P, rate) = LinearLoss(rate)(P)

"""
    QuadraticLoss(rate)

Quadratic mortality / loss functor.

!!! formulation
    ``l`` * ``P``²

    where:
    - ``P`` = plankton concentration
    - ``l`` = mortality (loss) rate

Quadratic mortality increases nonlinearly with biomass and is often interpreted to
represent viral processes and other density-dependent loss effects.
"""
struct QuadraticLoss{T}
    rate::T
end

@inline (q::QuadraticLoss)(P) = q.rate * P * P

"""
    quadratic_loss(P, rate)

Quadratic mortality (loss) rate.

!!! formulation
    ``l`` * ``P``²

    where:
    - ``P`` = plankton concentration
    - ``l`` = mortality (loss) rate

# Arguments
- `P`: plankton concentration
- `rate`: mortality (loss) rate
"""
@inline quadratic_loss(P, rate) = QuadraticLoss(rate)(P)

end # module
