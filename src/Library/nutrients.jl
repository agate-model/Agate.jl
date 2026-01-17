"""Building-block functors for nutrient limitation."""

module Nutrients

export MonodLimitation, LiebigMinimum

"""
    MonodLimitation(K)

Nutrient limitation functor `R ↦ R / (K + R)`.

`K` is the half-saturation constant.
"""
struct MonodLimitation{T}
    K::T
end

@inline (m::MonodLimitation)(R) = R / (m.K + R)

"""
    LiebigMinimum()

Return the minimum of its inputs.
"""
struct LiebigMinimum end

@inline (l::LiebigMinimum)(a, b) = ifelse(a < b, a, b)

@inline function (l::LiebigMinimum)(a, b, c, rest...)
    return l(l(a, b), c, rest...)
end

@inline function (l::LiebigMinimum)(values::NTuple{N,T}) where {N,T}
    m = values[1]
    @inbounds for i in 2:N
        m = l(m, values[i])
    end
    return m
end

end # module
