"""Building-block functors for nutrient limitation."""

module Nutrients

export monod_limitation, liebig_minimum

"""
    MonodLimitation(K)

Monod (Michaelis–Menten) nutrient limitation functor.

!!! formulation
    ``R`` / (``K`` + ``R``)

    where:
    - ``R`` = nutrient concentration (e.g. N, P, Si)
    - ``K`` = half-saturation constant

!!! tip
    This functional form is sometimes also used for predation (≈ Holling type II).
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

"""
    monod_limitation(R, K)

Monod (Michaelis–Menten) nutrient limitation.

!!! formulation
    ``R`` / (``K`` + ``R``)

    where:
    - ``R`` = nutrient concentration
    - ``K`` = half-saturation constant

# Arguments
- `R`: nutrient concentration
- `K`: nutrient half-saturation constant

!!! tip
    This is a thin, inlined alias around `MonodLimitation(K)(R)`.
"""
@inline monod_limitation(R, K) = MonodLimitation(K)(R)

"""
    LiebigMinimum()

Liebig's law of the minimum: return the minimum of nutrient limitation factors.

!!! formulation
    minimum(nutrient_limitations)

# Arguments
- `nutrient_limitations`: limitation factors (e.g. γᴺ, γᴾ, γˢⁱ) provided as positional
  arguments or as an `NTuple`.
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

"""
    liebig_minimum(a, b, rest...)
    liebig_minimum(values::NTuple)

Return the minimum value among the given limitation factors.

!!! formulation
    minimum(nutrient_limitations)

# Arguments
- `a, b, rest...`: limitation factors
- `values`: an `NTuple` of limitation factors

This is an explicit alias around `LiebigMinimum()` for clearer model code.
"""
@inline liebig_minimum(a, b) = LiebigMinimum()(a, b)

@inline liebig_minimum(a, b, c, rest...) = LiebigMinimum()(a, b, c, rest...)

@inline liebig_minimum(values::NTuple{N,T}) where {N,T} = LiebigMinimum()(values)

end # module
