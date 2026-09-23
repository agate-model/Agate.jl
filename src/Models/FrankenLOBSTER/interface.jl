using Adapt: adapt
import Adapt: adapt_structure

import OceanBioME: chlorophyll
import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

using OceanBioME.Models.NutrientsPlanktonDetritusModels: NutrientsPlanktonDetritus
import OceanBioME.Models.NutrientsPlanktonDetritusModels:
    dissolved_waste, inorganic_waste, nutrient_uptake, solid_waste, chlorophyll_ratio
import OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels: grazing

"""OceanBioME plankton component backed by one compiled Agate FrankenLOBSTER runtime."""
struct FrankenLOBSTERPlankton{
    Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers,PhytoplanktonTracers,ChlorophyllRatio
}
    runtime::Runtime
    chlorophyll_ratio::ChlorophyllRatio
end

function FrankenLOBSTERPlankton(
    runtime, owned::Tuple, exchange::Tuple=();
    phytoplankton_tracers=runtime.metadata.pft_entities.P,
    chlorophyll_ratio=1.31,
)
    owned_type = mapreduce(name -> typeof(Val(name)), (A, B) -> Union{A,B}, owned)
    return FrankenLOBSTERPlankton{
        typeof(runtime),owned,owned_type,exchange,phytoplankton_tracers,typeof(chlorophyll_ratio)
    }(runtime, chlorophyll_ratio)
end

@inline required_biogeochemical_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers}
) where {Runtime,OwnedTracers} = OwnedTracers

@inline required_biogeochemical_auxiliary_fields(
    ::FrankenLOBSTERPlankton{Runtime}
) where {Runtime} = required_biogeochemical_auxiliary_fields(Runtime)

@inline biogeochemical_drift_velocity(plankton::FrankenLOBSTERPlankton, tracer::Val) =
    biogeochemical_drift_velocity(plankton.runtime, tracer)

@inline chlorophyll_ratio(plankton::FrankenLOBSTERPlankton) = plankton.chlorophyll_ratio

@inline function chlorophyll(
    plankton::FrankenLOBSTERPlankton{R,O,T,E,P}, model
) where {R,O,T,E,P}
    biomass = mapreduce(name -> getproperty(model.tracers, name), +, P)
    return plankton.chlorophyll_ratio * biomass
end

@inline function adapt_structure(
    to, plankton::FrankenLOBSTERPlankton{R,O,T,E,P}
) where {R,O,T,E,P}
    return FrankenLOBSTERPlankton(
        adapt(to, plankton.runtime), O, E;
        phytoplankton_tracers=P,
        chlorophyll_ratio=adapt(to, plankton.chlorophyll_ratio),
    )
end

# OceanBioME fields -> Agate's statically ordered positional state.
@inline function _runtime_tracer_value(
    ::Val{Tracer},
    plankton::FrankenLOBSTERPlankton{R,O,T,E},
    i, j, k, fields,
) where {Tracer,R,O,T,E}
    Tracer in E && return zero(@inbounds getproperty(fields, first(O))[i, j, k])
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
    plankton::FrankenLOBSTERPlankton, tracer::Val, i, j, k, t, fields, auxiliary_fields
)
    tracer_values = _runtime_tracer_values(plankton, i, j, k, fields)
    auxiliary_values = _runtime_auxiliary_values(plankton, i, j, k, auxiliary_fields)
    x = zero(t)
    return plankton.runtime(tracer, x, x, x, t, tracer_values..., auxiliary_values...)
end

@inline _exchange_tendency(plankton, tracer, i, j, k, grid, fields, auxiliary_fields) =
    _agate_tendency(plankton, tracer, i, j, k, zero(eltype(grid)), fields, auxiliary_fields)

# The NPD call overload is restricted to the realized Agate-owned living tracer union, so
# OceanBioME nutrient/detritus/carbon/oxygen tracers keep their native dispatch.
@inline function (
    bgc::NutrientsPlanktonDetritus{FT,NUT,PLA}
)(i, j, k, grid, tracer::OwnedTracerType, clock, fields, auxiliary_fields) where {
    FT,NUT,Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers,
    PLA<:FrankenLOBSTERPlankton{Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers}
}
    return _agate_tendency(
        bgc.plankton, tracer, i, j, k, clock.time, fields, auxiliary_fields
    )
end

@inline nutrient_uptake(
    i, j, k, grid, ::Val{:NO₃}, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = -_exchange_tendency(plankton, Val(:NO₃), i, j, k, grid, fields, auxiliary_fields)

# Ammonium remains OceanBioME state for regeneration/nitrification but is not a
# FrankenLOBSTER phytoplankton substrate in this release.
@inline nutrient_uptake(
    i, j, k, grid, ::Val{:NH₄}, ::FrankenLOBSTERPlankton,
    ::NutrientsPlanktonDetritus{FT}, fields, auxiliary_fields,
) where FT = zero(FT)

@inline nutrient_uptake(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = nutrient_uptake(
    i, j, k, grid, Val(:NO₃), plankton, bgc, fields, auxiliary_fields
)

@inline solid_waste(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = _exchange_tendency(plankton, Val(:solid_waste), i, j, k, grid, fields, auxiliary_fields)

@inline dissolved_waste(
    i, j, k, grid, ::FrankenLOBSTERPlankton,
    ::NutrientsPlanktonDetritus{FT}, fields, auxiliary_fields,
) where FT = zero(FT)

@inline inorganic_waste(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = _exchange_tendency(
    plankton, Val(:inorganic_waste), i, j, k, grid, fields, auxiliary_fields
)

# DissolvedParticulate uses `grazing` for biological removal from organic-matter pools.
@inline grazing(
    i, j, k, grid, ::Val{:DOM}, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = -_exchange_tendency(plankton, Val(:DOM), i, j, k, grid, fields, auxiliary_fields)
