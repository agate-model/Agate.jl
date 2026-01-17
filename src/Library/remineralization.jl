"""Remineralization functors."""

module Remineralization

export LinearRemineralization

"""
    LinearRemineralization(rate)

Linear remineralization `D ↦ rate * D`.
"""
struct LinearRemineralization{T}
    rate::T
end

@inline (r::LinearRemineralization)(D) = r.rate * D

end # module
