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

"""
    QuadraticLoss(rate)

Quadratic loss `P ↦ rate * P^2`.
"""
struct QuadraticLoss{T}
    rate::T
end

@inline (q::QuadraticLoss)(P) = q.rate * P * P

"""
    linear_loss(P, rate)

Convenience wrapper for `LinearLoss(rate)(P)`.
"""
@inline linear_loss(P, rate) = LinearLoss(rate)(P)

"""
    quadratic_loss(P, rate)

Convenience wrapper for `QuadraticLoss(rate)(P)`.
"""
@inline quadratic_loss(P, rate) = QuadraticLoss(rate)(P)

end # module
