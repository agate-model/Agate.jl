"""Introspection helpers for constructed Agate biogeochemistry instances.

Tools to discover what a constructed model expects, without digging into the
generated types.
"""

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

@inline function _parameter_payload(params)
    return Base.hasproperty(params, :data) ? getproperty(params, :data) : params
end

"""    tracer_names(bgc) -> Vector{Symbol}

Return the ordered tracer symbols required by `bgc`.

The ordering matches Oceananigans / OceanBioME state-vector conventions.
"""
function tracer_names(bgc)::Vector{Symbol}
    return collect(required_biogeochemical_tracers(bgc))
end

"""    parameter_names(bgc) -> Vector{Symbol}

Return the parameter keys available on `bgc.parameters`.

Agate resolves only the parameters required by the constructed equations, so
this list is typically a *minimal* runtime set.
"""
function parameter_names(bgc)::Vector{Symbol}
    params = getproperty(bgc, :parameters)
    payload = _parameter_payload(params)
    names = collect(propertynames(payload))
    filter!(n -> n !== :data, names)
    return names
end

"""    required_parameters(bgc) -> Vector{Symbol}

Alias for [`parameter_names`](@ref).
"""
required_parameters(bgc)::Vector{Symbol} = parameter_names(bgc)
