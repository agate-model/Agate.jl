"""Building-block functors for nutrient limitation."""

module Nutrients

export monod_limitation, liebig_minimum, frank_minimum

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

@inline function (l::LiebigMinimum)(a, b)
    m = ifelse(a < b, a, b)
    return ifelse(isnan(a) | isnan(b), a + b, m)
end

@inline function (l::LiebigMinimum)(a, b, c, rest...)
    return l(l(a, b), c, rest...)
end

@inline function (l::LiebigMinimum)(values::Tuple{Vararg{Any,N}}) where N
    m = values[1]
    @inbounds for i in 2:N
        m = l(m, values[i])
    end
    return m
end

raw"""
    FrankMinimum()
    FrankMinimum(sharpness)

Differentiable Frank-family approximation to Liebig's minimum on limitation
factors in ``[0, 1]``.

!!! formulation
    For two limitation factors ``a`` and ``b``, let ``q = \exp(-s)`` where
    ``s`` is `sharpness`. The Frank minimum is

    ```math
    F_s(a, b) = \log_q\left(1 + \frac{(q^a - 1)(q^b - 1)}{q - 1}\right).
    ```

    Positive `sharpness` values give the minimum-like branch of the Frank
    family, with larger values approaching `min(a, b)`. The implementation uses
    an equivalent shifted form for numerical stability.

    ``1`` is the neutral element, so ``F_s(a, 1) = a``. ``0`` is absorbing, so
    ``F_s(a, 0) = 0``. For more than two factors the associative binary operator
    is applied successively.

Finite `sharpness` provides a smooth transition through nutrient co-limitation
for automatic differentiation; `LiebigMinimum` remains the exact hard minimum.
The default `sharpness = 50` keeps the transition localized while retaining a
smooth derivative crossover. Finite sharpness can underestimate the hard minimum
when several small limitation factors are similar; increasing `sharpness` reduces
that discrepancy while narrowing the smooth transition.

!!! warning "Domain"
    `FrankMinimum` is defined for normalized limitation factors in ``[0, 1]``.
    Inputs outside this interval can violate the minimum-like bounds or produce
    non-finite values.
"""
struct FrankMinimum{S}
    sharpness::S
end

const DEFAULT_FRANK_SHARPNESS = 50

@inline FrankMinimum() = FrankMinimum(DEFAULT_FRANK_SHARPNESS)

@inline function (f::FrankMinimum)(a, b)
    s = oftype(one(a + b), f.sharpness)
    a_is_min = a < b
    m = ifelse(a_is_min, a, b)
    M = ifelse(a_is_min, b, a)
    one_m = one(m)

    numerator = one_m + exp(-s * (M - m)) - exp(-s * M) - exp(-s * (one_m - m))
    denominator = one_m - exp(-s)
    result = m - log(numerator / denominator) / s

    return ifelse(isnan(a) | isnan(b), a + b, result)
end

@inline function (f::FrankMinimum)(a, b, c, rest...)
    return f((a, b, c, rest...))
end

@inline function (f::FrankMinimum)(values::Tuple{Vararg{Any,N}}) where N
    result = values[1]
    @inbounds for i in 2:N
        result = f(result, values[i])
    end
    return result
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

@inline liebig_minimum(values::Tuple{Vararg{Any,N}}) where N = LiebigMinimum()(values)

"""
    frank_minimum(a, b, rest...; sharpness = 50)
    frank_minimum(values::NTuple; sharpness = 50)

Return the differentiable Frank-family minimum of normalized limitation
factors in `[0, 1]`. Positive `sharpness` values increasingly approximate
`liebig_minimum`. Finite sharpness smooths co-limitation but can underestimate
the hard minimum when several small limitation factors are similar.
"""
@inline frank_minimum(a, b; sharpness=DEFAULT_FRANK_SHARPNESS) =
    FrankMinimum(sharpness)(a, b)

@inline frank_minimum(a, b, c, rest...; sharpness=DEFAULT_FRANK_SHARPNESS) =
    FrankMinimum(sharpness)(a, b, c, rest...)

@inline frank_minimum(
    values::Tuple{Vararg{Any,N}}; sharpness=DEFAULT_FRANK_SHARPNESS
) where N = FrankMinimum(sharpness)(values)

end # module
