module Remineralization

using ..Equations: AExpr, bgc_param

export remineralization_idealized
export remineralization_flux

"""
    remineralization_idealized(D, remineralization_rate)

Idealized remineralization of detritus into dissolved nutrients.

!!! formulation
    r * D

    where:
    - D = detritus concentration
    - r = remineralization rate

# Arguments
- `D`: detritus concentration
- `remineralization_rate`: remineralization rate
"""
@inline function remineralization_idealized(D, remineralization_rate)
    return remineralization_rate * D
end

"""\
    remineralization_flux(pool_sym::Symbol, rate_key::Symbol) -> AExpr

Construction-time symbolic helper for `rate_key * pool_sym`.

- `pool_sym` is a tracer symbol like `:DOC`.
- `rate_key` is a biogeochemical scalar key like `:DOC_remineralization`.

The scalar is recorded as an equation requirement and validated by the constructor.
"""
function remineralization_flux(pool_sym::Symbol, rate_key::Symbol)
    return bgc_param(rate_key) * pool_sym
end

end # module
