using Adapt: adapt
import Adapt: adapt_structure

import OceanBioME: chlorophyll
import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

using OceanBioME.Models.NutrientsPlanktonDetritusModels: NutrientsPlanktonDetritus
import OceanBioME.Models.NutrientsPlanktonDetritusModels:
    carbon_ratio,
    dissolved_waste,
    inorganic_waste,
    nutrient_uptake,
    solid_waste,
    chlorophyll_ratio
import OceanBioME.Models.NutrientsPlanktonDetritusModels.InorganicCarbonModels:
    biological_calcium_carbonate_dissolution,
    biological_calcium_carbonate_precipitation,
    particulate_calcium_carbonate_production
import OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels: grazing

"""OceanBioME plankton component backed by one compiled Agate FrankenLOBSTER runtime."""
struct FrankenLOBSTERPlankton{
    Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers,PhytoplanktonTracers,
    ProcessDiagnostics,T,
}
    runtime::Runtime
    process_diagnostics::ProcessDiagnostics
    chlorophyll_ratio::T
    carbon_ratio::T
    calcium_carbonate_rain_ratio::T
    zooplankton_calcium_carbonate_dissolution::T
end

function FrankenLOBSTERPlankton(
    runtime, owned::Tuple, exchange::Tuple;
    phytoplankton_tracers, process_diagnostics, chlorophyll_ratio, carbon_ratio,
    calcium_carbonate_rain_ratio, zooplankton_calcium_carbonate_dissolution,
)
    owned_type = mapreduce(name -> typeof(Val(name)), (A, B) -> Union{A,B}, owned)
    T = typeof(chlorophyll_ratio)
    return FrankenLOBSTERPlankton{
        typeof(runtime),owned,owned_type,exchange,phytoplankton_tracers,
        typeof(process_diagnostics),T,
    }(
        runtime, process_diagnostics, chlorophyll_ratio, carbon_ratio,
        calcium_carbonate_rain_ratio, zooplankton_calcium_carbonate_dissolution,
    )
end

@inline required_biogeochemical_tracers(
    ::FrankenLOBSTERPlankton{Runtime,OwnedTracers}
) where {Runtime,OwnedTracers} = OwnedTracers

@inline required_biogeochemical_tracers(
    npd::NutrientsPlanktonDetritus{FT,NUT,PLA}
) where {FT,NUT,PLA<:FrankenLOBSTERPlankton} = (
    required_biogeochemical_tracers(npd.nutrients)...,
    required_biogeochemical_tracers(npd.plankton)...,
    required_biogeochemical_tracers(npd.detritus)...,
    required_biogeochemical_tracers(npd.inorganic_carbon)...,
    required_biogeochemical_tracers(npd.oxygen)...,
    :T,
)

@inline required_biogeochemical_auxiliary_fields(
    ::FrankenLOBSTERPlankton{Runtime}
) where {Runtime} = required_biogeochemical_auxiliary_fields(Runtime)

@inline biogeochemical_drift_velocity(plankton::FrankenLOBSTERPlankton, tracer::Val) =
    biogeochemical_drift_velocity(plankton.runtime, tracer)

@inline chlorophyll_ratio(plankton::FrankenLOBSTERPlankton) = plankton.chlorophyll_ratio
@inline carbon_ratio(
    plankton::FrankenLOBSTERPlankton, ::NutrientsPlanktonDetritus{FT}
) where FT = convert(FT, plankton.carbon_ratio)

# P-specific calcite is supplied through the ExplicitCalciumCarbonate hooks below;
# CarbonateSystem therefore keeps OceanBioME's default implicit rain ratio.
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
        process_diagnostics=adapt(to, plankton.process_diagnostics),
        chlorophyll_ratio=adapt(to, plankton.chlorophyll_ratio),
        carbon_ratio=adapt(to, plankton.carbon_ratio),
        calcium_carbonate_rain_ratio=adapt(to, plankton.calcium_carbonate_rain_ratio),
        zooplankton_calcium_carbonate_dissolution=adapt(
            to, plankton.zooplankton_calcium_carbonate_dissolution
        ),
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
@inline (bgc::NutrientsPlanktonDetritus{<:Any,<:Any,PLA})(
    i, j, k, grid, tracer::OwnedTracerType, clock, fields, auxiliary_fields
) where {
    Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers,
    PLA<:FrankenLOBSTERPlankton{Runtime,OwnedTracers,OwnedTracerType,ExchangeTracers},
} = _agate_tendency(
    bgc.plankton, tracer, i, j, k, clock.time, fields, auxiliary_fields
)

@inline nutrient_uptake(
    i, j, k, grid, nutrient::Union{Val{:NO₃},Val{:NH₄}},
    plankton::FrankenLOBSTERPlankton, bgc::NutrientsPlanktonDetritus,
    fields, auxiliary_fields,
) = -_exchange_tendency(plankton, nutrient, i, j, k, grid, fields, auxiliary_fields)

@inline nutrient_uptake(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = nutrient_uptake(
    i, j, k, grid, Val(:NO₃), plankton, bgc, fields, auxiliary_fields
) + nutrient_uptake(
    i, j, k, grid, Val(:NH₄), plankton, bgc, fields, auxiliary_fields
)

@inline solid_waste(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = _exchange_tendency(plankton, Val(:solid_waste), i, j, k, grid, fields, auxiliary_fields)

@inline dissolved_waste(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = _exchange_tendency(
    plankton, Val(:dissolved_waste), i, j, k, grid, fields, auxiliary_fields
)

@inline inorganic_waste(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = _exchange_tendency(
    plankton, Val(:inorganic_waste), i, j, k, grid, fields, auxiliary_fields
)

@inline function _process_tendency(
    plankton::FrankenLOBSTERPlankton, ::Val{Process}, ::Val{Tracer},
    i, j, k, t, fields, auxiliary_fields,
) where {Process,Tracer}
    equations = getproperty(plankton.process_diagnostics, Process)
    hasfield(typeof(equations), Tracer) || return zero(t)
    equation = getfield(equations, Tracer)
    tracer_values = _runtime_tracer_values(plankton, i, j, k, fields)
    auxiliary_values = _runtime_auxiliary_values(plankton, i, j, k, auxiliary_fields)
    x = zero(t)
    return equation(plankton.runtime, x, x, x, t, tracer_values..., auxiliary_values...)
end

@inline function _phytoplankton_process_tendency(
    plankton::FrankenLOBSTERPlankton{R,O,T,E,P}, process::Val,
    i, j, k, grid, fields, auxiliary_fields,
) where {R,O,T,E,P}
    t = zero(eltype(grid))
    return mapreduce(
        tracer -> _process_tendency(
            plankton, process, Val(tracer), i, j, k, t, fields, auxiliary_fields
        ),
        +,
        P,
    )
end

@inline function biological_calcium_carbonate_precipitation(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
)
    retained_growth =
        _phytoplankton_process_tendency(
            plankton, Val(:nitrate_growth_P), i, j, k, grid, fields, auxiliary_fields
        ) +
        _phytoplankton_process_tendency(
            plankton, Val(:ammonia_growth_P), i, j, k, grid, fields, auxiliary_fields
        )
    return plankton.calcium_carbonate_rain_ratio * plankton.carbon_ratio * retained_growth
end

@inline function _phytoplankton_loss(
    plankton, process, i, j, k, grid, fields, auxiliary_fields
)
    return -_phytoplankton_process_tendency(
        plankton, process, i, j, k, grid, fields, auxiliary_fields
    )
end

@inline function particulate_calcium_carbonate_production(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
)
    grazing_loss = _phytoplankton_loss(
        plankton, Val(:grazing_Z_on_living), i, j, k, grid, fields, auxiliary_fields
    )
    mortality_loss = _phytoplankton_loss(
        plankton, Val(:mortality_P), i, j, k, grid, fields, auxiliary_fields
    )
    particulate_grazing =
        (one(plankton.zooplankton_calcium_carbonate_dissolution) -
         plankton.zooplankton_calcium_carbonate_dissolution) * grazing_loss
    return plankton.calcium_carbonate_rain_ratio * plankton.carbon_ratio *
           (particulate_grazing + mortality_loss)
end

@inline function biological_calcium_carbonate_dissolution(
    i, j, k, grid, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
)
    grazing_loss = _phytoplankton_loss(
        plankton, Val(:grazing_Z_on_living), i, j, k, grid, fields, auxiliary_fields
    )
    return plankton.calcium_carbonate_rain_ratio * plankton.carbon_ratio *
           plankton.zooplankton_calcium_carbonate_dissolution * grazing_loss
end

# DissolvedParticulate uses `grazing` for biological removal from organic-matter pools.
@inline grazing(
    i, j, k, grid, ::Val{:DOM}, plankton::FrankenLOBSTERPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) = -_exchange_tendency(plankton, Val(:DOM), i, j, k, grid, fields, auxiliary_fields)
