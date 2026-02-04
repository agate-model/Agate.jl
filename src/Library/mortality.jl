"""Mortality and loss functors."""

module Mortality

export linear_loss, quadratic_loss

"""
    LinearLoss(rate)

Linear loss `P ↦ rate * P`.
"""
struct LinearLoss{T}
    rate::T
end

@inline (l::LinearLoss)(P) = l.rate * P

"""Apply a linear loss rate to a state variable."""
@inline linear_loss(P, rate) = LinearLoss(rate)(P)

"""
    QuadraticLoss(rate)

Quadratic loss `P ↦ rate * P^2`.
"""
struct QuadraticLoss{T}
    rate::T
end

@inline (q::QuadraticLoss)(P) = q.rate * P * P

"""Apply a quadratic loss rate to a state variable."""
@inline quadratic_loss(P, rate) = QuadraticLoss(rate)(P)

end # module
