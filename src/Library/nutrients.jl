"""Nutrient-limitation and quota-regulation kernels."""

module Nutrients

export monod_limitation, liebig_minimum, frank_tnorm
export normalized_droop_limitation, quota_uptake_regulation

"""
    monod_limitation(R, K)

Return Monod (Michaelis-Menten) nutrient limitation ``R / (K + R)``.
The indeterminate `R == K == 0` case returns zero.
"""
@inline function monod_limitation(R, K)
    K == zero(K) && R == zero(R) && return zero(R)
    return R / (K + R)
end

"""
    liebig_minimum(a, b, rest...)
    liebig_minimum(values::NTuple)

Return the minimum of the supplied limitation factors while preserving NaN propagation.
"""
@inline function liebig_minimum(a, b)
    minimum_value = ifelse(a < b, a, b)
    return ifelse(isnan(a) | isnan(b), a + b, minimum_value)
end

@inline liebig_minimum(a, b, c, rest...) =
    liebig_minimum(liebig_minimum(a, b), c, rest...)

@inline function liebig_minimum(values::Tuple{Vararg{Any,N}}) where N
    result = values[1]
    @inbounds for i in 2:N
        result = liebig_minimum(result, values[i])
    end
    return result
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
