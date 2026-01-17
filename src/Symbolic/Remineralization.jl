module Remineralization

export remineralization_flux

"""
    remineralization_flux(PV, pool_sym, rate_key)

Construction-time symbolic helper for `PV[rate_key] * pool_sym`.
"""
@inline function remineralization_flux(PV, pool_sym::Symbol, rate_key::Symbol)
    rate = getproperty(PV, rate_key)
    return rate * pool_sym
end

end # module
