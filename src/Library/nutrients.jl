"""Building-block functors for nutrient limitation."""

module Nutrients

export monod_limitation, liebig_minimum, frank_tnorm
export normalized_droop_limitation, quota_uptake_regulation

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

const DEFAULT_FRANK_SHARPNESS = 50

@inline function _frank_tnorm_pair(a, b, sharpness)
    s = oftype(one(a + b), sharpness)
    a_is_min = a < b
    m = ifelse(a_is_min, a, b)
    M = ifelse(a_is_min, b, a)
    one_m = one(m)

    numerator = one_m + exp(-s * (M - m)) - exp(-s * M) - exp(-s * (one_m - m))
    denominator = one_m - exp(-s)
    result = m - log(numerator / denominator) / s

    return ifelse(isnan(a) | isnan(b), a + b, result)
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

`sharpness` is supplied explicitly as a numerical parameter; formulation
identity is represented separately by the process-authoring layer.
"""
@inline frank_tnorm(a, b; sharpness=DEFAULT_FRANK_SHARPNESS) =
    _frank_tnorm_pair(a, b, sharpness)

@inline frank_tnorm(a, b, c, rest...; sharpness=DEFAULT_FRANK_SHARPNESS) =
    frank_tnorm((a, b, c, rest...); sharpness)

@inline function frank_tnorm(
    values::Tuple{Vararg{Any,N}}; sharpness=DEFAULT_FRANK_SHARPNESS
) where N
    result = values[1]
    @inbounds for i in 2:N
        result = _frank_tnorm_pair(result, values[i], sharpness)
    end
    return result
end

"""
    normalized_droop_limitation(internal, reference, minimum_quota, maximum_quota)

Return a normalized Droop growth limitation from an internal nutrient inventory and
reference biomass inventory. The implied quota is `internal / reference`; limitation is
zero at or below `minimum_quota` and one at or above `maximum_quota`. Between those
bounds, the classical Droop response `1 - minimum_quota / quota` is normalized so that
`maximum_quota` maps exactly to one. Zero or negative reference biomass returns zero.
"""
@inline function normalized_droop_limitation(
    internal, reference, minimum_quota, maximum_quota
)
    reference > zero(reference) || return zero(internal + reference)
    quota = internal / reference
    quota > minimum_quota || return zero(quota)
    quota >= maximum_quota && return one(quota)
    response = (one(quota) - minimum_quota / quota) /
               (one(quota) - minimum_quota / maximum_quota)
    return min(one(response), max(zero(response), response))
end

"""
    quota_uptake_regulation(internal, reference, minimum_quota, maximum_quota, hill)

Return the bounded cellular-capacity factor used to regulate external nutrient uptake.
The response is one at or below `minimum_quota`, declines with normalized quota according
to `hill`, and is zero at or above `maximum_quota`. Zero or negative reference biomass
returns zero so an absent population has no uptake capacity.
"""
@inline function quota_uptake_regulation(
    internal, reference, minimum_quota, maximum_quota, hill
)
    reference > zero(reference) || return zero(internal + reference)
    quota = internal / reference
    quota <= minimum_quota && return one(quota)
    quota >= maximum_quota && return zero(quota)
    normalized = (quota - minimum_quota) / (maximum_quota - minimum_quota)
    response = one(normalized) - normalized^hill
    return min(one(response), max(zero(response), response))
end

end # module
