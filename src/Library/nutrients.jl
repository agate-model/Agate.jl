"""Building-block functors for nutrient limitation."""

module Nutrients

export monod_limitation, liebig_minimum, frank_tnorm

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

"""
    FrankTNorm()
    FrankTNorm(sharpness)

Differentiable Frank t-norm approximation to Liebig's minimum on limitation
factors in `[0, 1]`.

!!! formulation
    For two limitation factors `a` and `b`, let `q = exp(-s)`, where `s` is
    `sharpness`. The Frank t-norm is

    ```text
    q = exp(-s)
    F(a, b) = log(1 + ((q^a - 1) * (q^b - 1)) / (q - 1)) / log(q)
    ```

    Positive `sharpness` values give the minimum-like branch of the Frank
    family, with larger values approaching `min(a, b)`. The implementation uses
    an equivalent shifted form for numerical stability.

    `1` is the neutral element, so `F(a, 1) = a`. `0` is absorbing, so
    `F(a, 0) = 0`. For more than two factors the associative binary operator is
    applied successively.

Finite `sharpness` provides a smooth transition through nutrient co-limitation
for automatic differentiation; `LiebigMinimum` remains the exact hard minimum.
The default `sharpness = 50` keeps the transition localized while retaining a
smooth derivative crossover. Finite sharpness can underestimate the hard minimum
when several small limitation factors are similar; increasing `sharpness` reduces
that discrepancy while narrowing the smooth transition.

!!! warning "Domain"
    `FrankTNorm` is defined for normalized limitation factors in `[0, 1]`.
    Inputs outside this interval can violate the minimum-like bounds or produce
    non-finite values.
"""
struct FrankTNorm{S}
    sharpness::S
end

const DEFAULT_FRANK_SHARPNESS = 50

@inline FrankTNorm() = FrankTNorm(DEFAULT_FRANK_SHARPNESS)

@inline function (f::FrankTNorm)(a, b)
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

@inline function (f::FrankTNorm)(a, b, c, rest...)
    return f((a, b, c, rest...))
end

@inline function (f::FrankTNorm)(values::Tuple{Vararg{Any,N}}) where N
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
    frank_tnorm(a, b, rest...; sharpness = 50)
    frank_tnorm(values::NTuple; sharpness = 50)

Return the differentiable Frank t-norm approximation to Liebig's minimum for
normalized limitation factors in `[0, 1]`.

For two limitation factors `a` and `b`, let `q = exp(-s)`, where `s` is
`sharpness`. The Frank t-norm is

```text
q = exp(-s)
F(a, b) = log(1 + ((q^a - 1) * (q^b - 1)) / (q - 1)) / log(q)
```

Positive `sharpness` values give the minimum-like branch of the Frank family,
with larger values approaching `liebig_minimum(a, b)`. The implementation uses
an equivalent shifted form for numerical stability. For more than two factors,
the associative binary operator is applied successively.

`1` is the neutral element and `0` is absorbing. Finite sharpness provides a
smooth derivative transition through nutrient co-limitation for automatic
differentiation. The default `sharpness = 50` keeps that transition localized
while retaining a smooth crossover. Finite sharpness can underestimate the hard
minimum when several small limitation factors are similar; increasing
`sharpness` reduces that discrepancy while narrowing the smooth transition.

!!! warning "Domain"
    The Frank t-norm is defined for normalized limitation factors in `[0, 1]`.
    Inputs outside this interval can violate the minimum-like bounds or produce
    non-finite values.

The callable `FrankTNorm(sharpness)` functor provides the same operation for
lower-level tendency configuration.
"""
@inline frank_tnorm(a, b; sharpness=DEFAULT_FRANK_SHARPNESS) =
    FrankTNorm(sharpness)(a, b)

@inline frank_tnorm(a, b, c, rest...; sharpness=DEFAULT_FRANK_SHARPNESS) =
    FrankTNorm(sharpness)(a, b, c, rest...)

@inline frank_tnorm(
    values::Tuple{Vararg{Any,N}}; sharpness=DEFAULT_FRANK_SHARPNESS
) where N = FrankTNorm(sharpness)(values)

end # module
