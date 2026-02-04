"""Building-block functors for nutrient limitation."""

module Nutrients

export monod_limitation, liebig_minimum

"""
    MonodLimitation(K)

Nutrient limitation functor `R ↦ R / (K + R)`.

`K` is the half-saturation constant.
"""
struct MonodLimitation{T}
    K::T
end

@inline function (m::MonodLimitation)(R)
    K = m.K
    if K == zero(K) && R == zero(R)
        return zero(R)
    end
    return R / (K + R)
end

"""Apply Monod nutrient limitation `R/(K+R)`.

This is a thin, inlined alias around `MonodLimitation(K)(R)` for clearer model code.
"""
@inline monod_limitation(R, K) = MonodLimitation(K)(R)

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

"""Return the Liebig minimum across inputs.

This is an explicit alias around `LiebigMinimum()` for clearer model code.
"""
@inline liebig_minimum(a, b) = LiebigMinimum()(a, b)

@inline liebig_minimum(a, b, c, rest...) = LiebigMinimum()(a, b, c, rest...)

@inline liebig_minimum(values::NTuple{N,T}) where {N,T} = LiebigMinimum()(values)

end # module
