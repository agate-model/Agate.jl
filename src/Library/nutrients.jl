"""Building-block functors for nutrient limitation."""

module Nutrients

export monod_limitation, liebig_minimum, smooth_liebig_minimum

"""
    MonodLimitation(K)

Monod (Michaelis–Menten) nutrient limitation functor.

!!! formulation
    ``R`` / (``K`` + ``R``)

    where:
    - ``R`` = nutrient concentration (e.g. N, P, Si)
    - ``K`` = half-saturation constant

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
    - ``K`` = nutrient half-saturation constant

# Arguments
- `R`: nutrient concentration
- `K`: nutrient half-saturation constant

!!! tip
    This functional form is sometimes also used for predation (≈ Holling type II).
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

raw"""
    SmoothLiebigMinimum(sharpness)

Smooth approximation to Liebig's law of the minimum.

!!! formulation
    ```math
    -\frac{1}{s}\log\sum_i \exp(-s x_i)
    ```

    where ``s`` is `sharpness`. Larger `sharpness` values more closely
    approximate the hard minimum. The implementation uses a shifted log-sum-exp
    form for numerical stability.
"""
struct SmoothLiebigMinimum{S}
    sharpness::S
end

@inline SmoothLiebigMinimum() = SmoothLiebigMinimum(50.0)

@inline function (l::SmoothLiebigMinimum)(a, b)
    s = l.sharpness
    m = ifelse(a < b, a, b)
    return m - log(exp(-s * (a - m)) + exp(-s * (b - m))) / s
end

@inline function (l::SmoothLiebigMinimum)(a, b, c, rest...)
    return l((a, b, c, rest...))
end

@inline function (l::SmoothLiebigMinimum)(values::NTuple{N,T}) where {N,T}
    s = l.sharpness
    m = values[1]
    @inbounds for i in 2:N
        m = ifelse(values[i] < m, values[i], m)
    end

    total = zero(m)
    @inbounds for i in 1:N
        total += exp(-s * (values[i] - m))
    end

    return m - log(total) / s
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

"""
    smooth_liebig_minimum(a, b, rest...; sharpness = 50.0)
    smooth_liebig_minimum(values::NTuple; sharpness = 50.0)

Return a smooth approximation to the minimum value among the given limitation
factors. Larger `sharpness` values approach `liebig_minimum`.
"""
@inline smooth_liebig_minimum(a, b; sharpness=50.0) = SmoothLiebigMinimum(sharpness)(a, b)

@inline smooth_liebig_minimum(a, b, c, rest...; sharpness=50.0) = SmoothLiebigMinimum(sharpness)(a, b, c, rest...)

@inline smooth_liebig_minimum(values::NTuple{N,T}; sharpness=50.0) where {N,T} = SmoothLiebigMinimum(sharpness)(values)

end # module
