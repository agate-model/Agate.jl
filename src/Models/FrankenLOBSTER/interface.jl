using Adapt: adapt
import Adapt: adapt_structure

import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus
import OceanBioME.Models.NutrientsPlanktonDetritusModels:
    dissolved_waste,
    inorganic_waste,
    nutrient_uptake,
    solid_waste
import OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels: grazing

"""Internal OceanBioME plankton component backed by a compiled Agate realization.

`OwnedTracers` are living P/Z/B tracers registered by NPD and `ExchangeTracers` are
compiled accumulators reported through NPD hooks rather than registered as prognostic fields.
OceanBioME-owned external fields are the remaining tracer identities in the compiled runtime.
"""
struct FrankenLOBSTERPlankton{Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers}
    runtime::Runtime
end

function FrankenLOBSTERPlankton(runtime, owned::Tuple, exchange::Tuple=())
    owned_type = mapreduce(name -> typeof(Val(name)), (A, B) -> Union{A,B}, owned)
    return FrankenLOBSTERPlankton{typeof(runtime),owned,owned_type,exchange}(runtime)
end

@inline required_biogeochemical_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers}
) where {Runtime,OwnedTracers} = OwnedTracers

@inline required_biogeochemical_auxiliary_fields(
    ::FrankenLOBSTERPlankton{Runtime}
) where {Runtime} = required_biogeochemical_auxiliary_fields(Runtime)

@inline biogeochemical_drift_velocity(plankton::FrankenLOBSTERPlankton, tracer::Val) =
    biogeochemical_drift_velocity(plankton.runtime, tracer)

@inline function external_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers}
) where {Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers}
    return Tuple(
        tracer for tracer in required_biogeochemical_tracers(Runtime)
        if tracer ∉ OwnedTracers && tracer ∉ ExchangeTracers
    )
end

@inline exchange_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers}
) where {Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers} = ExchangeTracers

@inline function adapt_structure(
    to,
    plankton::FrankenLOBSTERPlankton{Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers},
) where {Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers}
    return FrankenLOBSTERPlankton(adapt(to, plankton.runtime), OwnedTracers, ExchangeTracers)
end

# Static OceanBioME field -> Agate positional-runtime bridge. AgateBGC already type-encodes
# its complete tracer and auxiliary identity tuples, so the adapter reuses those directly.
@inline function _exchange_zero(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers}, i, j, k, fields
) where {Runtime,OwnedTracers}
    return zero(@inbounds getproperty(fields, first(OwnedTracers))[i, j, k])
end

@inline function _runtime_tracer_value(
    ::Val{Tracer},
    plankton::FrankenLOBSTERPlankton{Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers},
    i,
    j,
    k,
    fields,
) where {Tracer,Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers}
    Tracer in ExchangeTracers && return _exchange_zero(plankton, i, j, k, fields)
    return @inbounds getproperty(fields, Tracer)[i, j, k]
end

@inline function _runtime_tracer_values(
    plankton::FrankenLOBSTERPlankton{Runtime}, i, j, k, fields
) where {Runtime}
    tracers = required_biogeochemical_tracers(Runtime)
    return ntuple(Val(length(tracers))) do n
        _runtime_tracer_value(Val(tracers[n]), plankton, i, j, k, fields)
    end
end

@inline function _runtime_auxiliary_values(
    ::FrankenLOBSTERPlankton{Runtime}, i, j, k, auxiliary_fields
) where {Runtime}
    auxiliaries = required_biogeochemical_auxiliary_fields(Runtime)
    return ntuple(Val(length(auxiliaries))) do n
        @inbounds getproperty(auxiliary_fields, auxiliaries[n])[i, j, k]
    end
end

@inline function _agate_tendency(
    plankton::FrankenLOBSTERPlankton,
    tracer::Val,
    i,
    j,
    k,
    clock,
    fields,
    auxiliary_fields,
)
    tracer_values = _runtime_tracer_values(plankton, i, j, k, fields)
    auxiliary_values = _runtime_auxiliary_values(plankton, i, j, k, auxiliary_fields)
    t = clock.time
    x = zero(t)
    return plankton.runtime(tracer, x, x, x, t, tracer_values..., auxiliary_values...)
end

@inline _zero_clock(grid) = (; time=zero(eltype(grid)))

# Restrict the NPD call overload to the concrete union of Agate-owned living tracer
# Val types encoded in the plankton wrapper. OceanBioME-owned nutrient, detritus,
# carbon, and oxygen tracers therefore keep their native NPD dispatch unchanged.
@inline function (
    bgc::NutrientsPlanktonDetritus{FT,NUT,PLA}
)(i, j, k, grid, tracer::OwnedTracerType, clock, fields, auxiliary_fields) where {
    FT,
    NUT,
    Runtime,
    OwnedTracers,
    OwnedTracerType,
    ExchangeTracers,
    PLA<:FrankenLOBSTERPlankton{
        Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers
    },
}
    return _agate_tendency(
        bgc.plankton, tracer, i, j, k, clock, fields, auxiliary_fields
    )
end

@inline nutrient_uptake(
    i,
    j,
    k,
    grid,
    tracer::Union{Val{:NO₃},Val{:NH₄}},
    plankton::FrankenLOBSTERPlankton,
    bgc,
    fields,
    auxiliary_fields,
) = -_agate_tendency(
    plankton, tracer, i, j, k, _zero_clock(grid), fields, auxiliary_fields
)

@inline solid_waste(
    i,
    j,
    k,
    grid,
    plankton::FrankenLOBSTERPlankton,
    bgc,
    fields,
    auxiliary_fields,
) = _agate_tendency(
    plankton,
    Val(:solid_waste),
    i,
    j,
    k,
    _zero_clock(grid),
    fields,
    auxiliary_fields,
)

@inline dissolved_waste(
    i,
    j,
    k,
    grid,
    ::FrankenLOBSTERPlankton,
    ::NutrientsPlanktonDetritus{FT},
    fields,
    auxiliary_fields,
) where FT = zero(FT)

@inline inorganic_waste(
    i,
    j,
    k,
    grid,
    plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus,
    fields,
    auxiliary_fields,
) = _agate_tendency(
    plankton,
    Val(:inorganic_waste),
    i,
    j,
    k,
    _zero_clock(grid),
    fields,
    auxiliary_fields,
)

# NPD's organic-matter components call `grazing` for biological removal. The hook is
# resource-generic on the Agate side; DOM is the first active FrankenLOBSTER substrate.
@inline grazing(
    i,
    j,
    k,
    grid,
    ::Val{:DOM},
    plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus{FT},
    fields,
    auxiliary_fields,
) where FT = -_agate_tendency(
    plankton, Val(:DOM), i, j, k, _zero_clock(grid), fields, auxiliary_fields
)
