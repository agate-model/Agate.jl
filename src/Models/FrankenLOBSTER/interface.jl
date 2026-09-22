using Adapt: adapt
import Adapt: adapt_structure

import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

using OceanBioME.Models.NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus
using OceanBioME.Models.NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients,
    NitrateAmmonia

import OceanBioME.Models.NutrientsPlanktonDetritusModels:
    dissolved_waste,
    inorganic_nitrogen_waste,
    inorganic_waste,
    nutrient_uptake,
    solid_waste

"""Internal OceanBioME plankton component backed by a compiled Agate realization.

`OwnedTracers` are living P/Z/B tracers registered by NPD and `ExchangeTracers` are
compiled accumulators reported through NPD hooks rather than registered as prognostic fields.
OceanBioME-owned external fields are the remaining tracer identities in the compiled runtime.
"""
struct FrankenLOBSTERPlankton{Runtime,OwnedTracers,ExchangeTracers}
    runtime::Runtime
end

FrankenLOBSTERPlankton(runtime, owned::Tuple, exchange::Tuple=()) =
    FrankenLOBSTERPlankton{typeof(runtime),owned,exchange}(runtime)

@inline required_biogeochemical_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers}
) where {Runtime,OwnedTracers} = OwnedTracers

@inline required_biogeochemical_auxiliary_fields(
    ::FrankenLOBSTERPlankton{Runtime}
) where {Runtime} = required_biogeochemical_auxiliary_fields(Runtime)

@inline biogeochemical_drift_velocity(plankton::FrankenLOBSTERPlankton, tracer::Val) =
    biogeochemical_drift_velocity(plankton.runtime, tracer)

@inline function external_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers,ExchangeTracers}
) where {Runtime,OwnedTracers,ExchangeTracers}
    return Tuple(
        tracer for tracer in required_biogeochemical_tracers(Runtime)
        if tracer ∉ OwnedTracers && tracer ∉ ExchangeTracers
    )
end

@inline exchange_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers,ExchangeTracers}
) where {Runtime,OwnedTracers,ExchangeTracers} = ExchangeTracers

@inline function adapt_structure(
    to,
    plankton::FrankenLOBSTERPlankton{Runtime,OwnedTracers,ExchangeTracers},
) where {Runtime,OwnedTracers,ExchangeTracers}
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
    plankton::FrankenLOBSTERPlankton{Runtime,OwnedTracers,ExchangeTracers},
    i,
    j,
    k,
    fields,
) where {Tracer,Runtime,OwnedTracers,ExchangeTracers}
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

# NPD plankton contract. The generic tendency method is used for arbitrary living-tracer
# names; exact nitrate/ammonium intersections below preserve OceanBioME nutrient dispatch.
@inline function (
    bgc::NutrientsPlanktonDetritus{FT,NUT,PLA}
)(i, j, k, grid, tracer::Val, clock, fields, auxiliary_fields) where {
    FT,NUT,PLA<:FrankenLOBSTERPlankton
}
    return _agate_tendency(bgc.plankton, tracer, i, j, k, clock, fields, auxiliary_fields)
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
    ::FrankenLOBSTERPlankton,
    ::NutrientsPlanktonDetritus{FT},
    fields,
    auxiliary_fields,
) where FT = zero(FT)

@inline function (
    bgc::NutrientsPlanktonDetritus{FT,NUT,PLA}
)(i, j, k, grid, tracer::Val{:NO₃}, clock, fields, auxiliary_fields) where {
    FT,
    NUT<:Nutrients{<:NitrateAmmonia},
    PLA<:FrankenLOBSTERPlankton,
}
    nitrification = @inbounds fields.NH₄[i, j, k] * bgc.nutrients.nitrogen.nitrification_rate
    return nitrification - nutrient_uptake(
        i, j, k, grid, tracer, bgc.plankton, bgc, fields, auxiliary_fields
    )
end

@inline function (
    bgc::NutrientsPlanktonDetritus{FT,NUT,PLA}
)(i, j, k, grid, tracer::Val{:NH₄}, clock, fields, auxiliary_fields) where {
    FT,
    NUT<:Nutrients{<:NitrateAmmonia},
    PLA<:FrankenLOBSTERPlankton,
}
    nitrification = @inbounds fields.NH₄[i, j, k] * bgc.nutrients.nitrogen.nitrification_rate
    regenerated =
        inorganic_nitrogen_waste(
            i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields
        ) +
        inorganic_nitrogen_waste(
            i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields
        )
    return regenerated - nutrient_uptake(
        i, j, k, grid, tracer, bgc.plankton, bgc, fields, auxiliary_fields
    ) - nitrification
end
