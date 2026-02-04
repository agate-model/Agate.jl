"""Remineralization functors."""

module Remineralization

export linear_remineralization

"""
    LinearRemineralization(rate)

Linear remineralization `D ↦ rate * D`.
"""
struct LinearRemineralization{T}
    rate::T
end

@inline (r::LinearRemineralization)(D) = r.rate * D

"""Apply a linear remineralization rate to detrital material."""
@inline linear_remineralization(D, rate) = LinearRemineralization(rate)(D)

end # module
