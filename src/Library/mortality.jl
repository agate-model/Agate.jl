"""Mortality and loss functors."""

module Mortality

export LinearLoss, QuadraticLoss

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

end # module
