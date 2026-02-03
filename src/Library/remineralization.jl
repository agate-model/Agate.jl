"""Remineralization functors."""

module Remineralization

export remin

"""
    LinearRemineralization(rate)

Linear remineralization `D ↦ rate * D`.
"""
struct LinearRemineralization{T}
    rate::T
end

@inline (r::LinearRemineralization)(D) = r.rate * D

"""
    remin(D, rate)

Convenience wrapper for `LinearRemineralization(rate)(D)`.
"""
@inline remin(D, rate) = LinearRemineralization(rate)(D)

end # module
