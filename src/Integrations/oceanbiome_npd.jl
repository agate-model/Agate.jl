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
    chlorophyll_ratio,
    dissolved_waste,
    inorganic_waste,
    nutrient_uptake,
    solid_waste
import OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels: grazing

"""OceanBioME NPD plankton component backed by one compiled Agate runtime.

`NPDPlankton` is an integration boundary, not a biological model. The compiled
Agate runtime owns the ecological equations; this wrapper maps their signed
tracer tendencies onto OceanBioME's existing NPD plankton hooks.
"""
struct NPDPlankton{
    Coupling,
    Runtime,
    OwnedTracers,
    OwnedTracerType,
    NutrientTracers,
    ExchangeTracers,
    ConsumedDetritus,
    Dependencies,
    PhytoplanktonTracers,
    Traits,
}
    runtime::Runtime
    traits::Traits
    coupling::Coupling
end

function _component_tracers(runtime, components::Tuple)
    metadata = runtime.metadata
    hasproperty(metadata, :component_tracers) || throw(
        ArgumentError("Agate runtime metadata does not expose component tracer identities."),
    )
    return Tuple(
        tracer
        for component in components
        for tracer in begin
            hasproperty(metadata.component_tracers, component) || throw(
                ArgumentError("Unknown Agate component :$component."),
            )
            getproperty(metadata.component_tracers, component)
        end
    )
end

function _validate_npd_traits(traits::NamedTuple)
    for name in (:carbon_ratio, :chlorophyll_ratio)
        hasproperty(traits, name) || throw(
            ArgumentError("NPDPlankton traits must define :$name."),
        )
        value = getproperty(traits, name)
        value isa Real && !(value isa Bool) && isfinite(value) && value >= 0 || throw(
            ArgumentError("NPDPlankton trait :$name must be finite and nonnegative."),
        )
    end
    traits.carbon_ratio > 0 || throw(ArgumentError("NPDPlankton carbon_ratio must be > 0."))
    return traits
end

"""
    NPDPlankton(runtime; owned_components, phytoplankton_components=(), nutrient_tracers=(),
                exchange_tracers=(solid=:solid_waste, dissolved=:dissolved_waste,
                                  inorganic=:inorganic_waste),
                consumed_detritus=(), dependencies=(), traits, coupling=nothing)

Wrap a compiled Agate runtime as an OceanBioME `NutrientsPlanktonDetritus` plankton component.
Component names are resolved once from Agate runtime metadata; all cell-level coupling is then
statically dispatched from the resulting tracer tuples.
"""
function NPDPlankton(
    runtime;
    owned_components::Tuple,
    phytoplankton_components::Tuple=(),
    nutrient_tracers::Tuple=(),
    exchange_tracers::NamedTuple=(
        solid=:solid_waste,
        dissolved=:dissolved_waste,
        inorganic=:inorganic_waste,
    ),
    consumed_detritus::Tuple=(),
    dependencies::Tuple=(),
    traits::NamedTuple,
    coupling=nothing,
)
    keys(exchange_tracers) == (:solid, :dissolved, :inorganic) || throw(
        ArgumentError("exchange_tracers must define (:solid, :dissolved, :inorganic)."),
    )
    all(name -> name in (:NO₃, :NH₄), nutrient_tracers) || throw(
        ArgumentError("NPDPlankton nutrient_tracers currently supports only :NO₃ and :NH₄."),
    )
    owned = _component_tracers(runtime, owned_components)
    isempty(owned) && throw(ArgumentError("NPDPlankton must own at least one tracer."))
    phytoplankton = _component_tracers(runtime, phytoplankton_components)
    isempty(phytoplankton) && throw(
        ArgumentError("NPDPlankton must identify at least one phytoplankton tracer."),
    )
    owned_type = mapreduce(name -> typeof(Val(name)), (A, B) -> Union{A,B}, owned)
    exchanges = Tuple(values(exchange_tracers))
    _validate_npd_traits(traits)

    return NPDPlankton{
        typeof(coupling),
        typeof(runtime),
        owned,
        owned_type,
        nutrient_tracers,
        exchanges,
        consumed_detritus,
        dependencies,
        phytoplankton,
        typeof(traits),
    }(runtime, traits, coupling)
end

@inline _owned_tracers(::NPDPlankton{C,R,O}) where {C,R,O} = O
@inline _nutrient_tracers(::NPDPlankton{C,R,O,OT,N}) where {C,R,O,OT,N} = N
@inline _exchange_tracers(::NPDPlankton{C,R,O,OT,N,E}) where {C,R,O,OT,N,E} = E
@inline _consumed_detritus(::NPDPlankton{C,R,O,OT,N,E,D}) where {C,R,O,OT,N,E,D} = D
@inline _dependencies(::NPDPlankton{C,R,O,OT,N,E,D,Deps}) where {C,R,O,OT,N,E,D,Deps} = Deps
@inline phytoplankton_tracers(::NPDPlankton{C,R,O,OT,N,E,D,Deps,P}) where {C,R,O,OT,N,E,D,Deps,P} = P

@inline required_biogeochemical_tracers(plankton::NPDPlankton) = _owned_tracers(plankton)
@inline required_biogeochemical_auxiliary_fields(plankton::NPDPlankton) =
    required_biogeochemical_auxiliary_fields(plankton.runtime)
@inline biogeochemical_drift_velocity(plankton::NPDPlankton, tracer::Val) =
    biogeochemical_drift_velocity(plankton.runtime, tracer)

@inline chlorophyll_ratio(plankton::NPDPlankton) = plankton.traits.chlorophyll_ratio
@inline carbon_ratio(plankton::NPDPlankton, ::NutrientsPlanktonDetritus{FT}) where FT =
    convert(FT, plankton.traits.carbon_ratio)
@inline chlorophyll(plankton::NPDPlankton, model) = plankton.traits.chlorophyll_ratio *
    mapreduce(name -> getproperty(model.tracers, name), +, phytoplankton_tracers(plankton))

@inline function adapt_structure(to, plankton::NPDPlankton{C,R,O,OT,N,E,D,Deps,P,T}) where {C,R,O,OT,N,E,D,Deps,P,T}
    runtime = adapt(to, plankton.runtime)
    traits = adapt(to, plankton.traits)
    coupling = adapt(to, plankton.coupling)
    return NPDPlankton{typeof(coupling),typeof(runtime),O,OT,N,E,D,Deps,P,typeof(traits)}(
        runtime, traits, coupling
    )
end

@inline function _append_unique(acc::Tuple, values::Tuple)
    isempty(values) && return acc
    head = first(values)
    next = head in acc ? acc : (acc..., head)
    return _append_unique(next, Base.tail(values))
end

@inline function required_biogeochemical_tracers(
    npd::NutrientsPlanktonDetritus{FT,NUT,PLA},
) where {FT,NUT,PLA<:NPDPlankton}
    tracers = (
        required_biogeochemical_tracers(npd.nutrients)...,
        required_biogeochemical_tracers(npd.plankton)...,
        required_biogeochemical_tracers(npd.detritus)...,
        required_biogeochemical_tracers(npd.inorganic_carbon)...,
        required_biogeochemical_tracers(npd.oxygen)...,
        _dependencies(npd.plankton)...,
    )
    return _append_unique((), tracers)
end

# OceanBioME fields -> Agate's statically ordered positional state.
@inline function _runtime_tracer_value(::Val{Tracer}, plankton::NPDPlankton, i, j, k, fields) where Tracer
    Tracer in _exchange_tracers(plankton) &&
        return zero(@inbounds getproperty(fields, first(_owned_tracers(plankton)))[i, j, k])
    return @inbounds getproperty(fields, Tracer)[i, j, k]
end

@inline function _runtime_tracer_values(plankton::NPDPlankton, i, j, k, fields)
    tracers = required_biogeochemical_tracers(plankton.runtime)
    return ntuple(Val(length(tracers))) do n
        _runtime_tracer_value(Val(tracers[n]), plankton, i, j, k, fields)
    end
end

@inline function _runtime_auxiliary_values(plankton::NPDPlankton, i, j, k, auxiliary_fields)
    auxiliaries = required_biogeochemical_auxiliary_fields(plankton.runtime)
    return ntuple(Val(length(auxiliaries))) do n
        @inbounds getproperty(auxiliary_fields, auxiliaries[n])[i, j, k]
    end
end

@inline function _agate_tendency(
    plankton::NPDPlankton, tracer::Val, i, j, k, t, fields, auxiliary_fields
)
    tracer_values = _runtime_tracer_values(plankton, i, j, k, fields)
    auxiliary_values = _runtime_auxiliary_values(plankton, i, j, k, auxiliary_fields)
    x = zero(t)
    return plankton.runtime(tracer, x, x, x, t, tracer_values..., auxiliary_values...)
end

@inline _exchange_tendency(plankton, tracer, i, j, k, grid, fields, auxiliary_fields) =
    _agate_tendency(plankton, tracer, i, j, k, zero(eltype(grid)), fields, auxiliary_fields)

# Restrict the NPD call overload to the Agate-owned living tracer union so OceanBioME
# nutrient/detritus/carbon/oxygen tracers keep their native dispatch.
@inline (bgc::NutrientsPlanktonDetritus{<:Any,<:Any,PLA})(
    i, j, k, grid, tracer::OwnedTracerType, clock, fields, auxiliary_fields
) where {
    C,Runtime,OwnedTracers,OwnedTracerType,N,E,D,Deps,P,T,
    PLA<:NPDPlankton{C,Runtime,OwnedTracers,OwnedTracerType,N,E,D,Deps,P,T},
} = _agate_tendency(
    bgc.plankton, tracer, i, j, k, clock.time, fields, auxiliary_fields
)

@inline function nutrient_uptake(
    i, j, k, grid, nutrient::Union{Val{:NO₃},Val{:NH₄}}, plankton::NPDPlankton,
    ::NutrientsPlanktonDetritus, fields, auxiliary_fields,
)
    name = nutrient isa Val{:NO₃} ? :NO₃ : :NH₄
    name in _nutrient_tracers(plankton) || return zero(eltype(grid))
    return -_exchange_tendency(plankton, nutrient, i, j, k, grid, fields, auxiliary_fields)
end

@inline _sum_nutrient_uptake(
    ::Tuple{}, i, j, k, grid, plankton, bgc, fields, auxiliary_fields,
) = zero(eltype(grid))

@inline function _sum_nutrient_uptake(
    nutrients::Tuple, i, j, k, grid, plankton, bgc, fields, auxiliary_fields,
)
    nutrient = first(nutrients)
    return nutrient_uptake(
        i, j, k, grid, Val(nutrient), plankton, bgc, fields, auxiliary_fields
    ) + _sum_nutrient_uptake(
        Base.tail(nutrients), i, j, k, grid, plankton, bgc, fields, auxiliary_fields
    )
end

@inline function nutrient_uptake(
    i, j, k, grid, plankton::NPDPlankton,
    bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
)
    return _sum_nutrient_uptake(
        _nutrient_tracers(plankton), i, j, k, grid, plankton, bgc, fields, auxiliary_fields
    )
end

@inline function _exchange_channel(
    plankton, channel, i, j, k, grid, fields, auxiliary_fields,
)
    channel === nothing && return zero(eltype(grid))
    return _exchange_tendency(
        plankton, Val(channel), i, j, k, grid, fields, auxiliary_fields
    )
end

for (hook, index) in ((:solid_waste, 1), (:dissolved_waste, 2), (:inorganic_waste, 3))
    @eval @inline function $hook(
        i, j, k, grid, plankton::NPDPlankton,
        ::NutrientsPlanktonDetritus, fields, auxiliary_fields,
    )
        return _exchange_channel(
            plankton, _exchange_tracers(plankton)[$index], i, j, k, grid, fields, auxiliary_fields
        )
    end
end

# DissolvedParticulate uses `grazing` for biological removal from organic-matter pools.
@inline function grazing(
    i, j, k, grid, ::Val{Tracer}, plankton::NPDPlankton,
    ::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) where Tracer
    Tracer in _consumed_detritus(plankton) || return zero(eltype(grid))
    return -_exchange_tendency(plankton, Val(Tracer), i, j, k, grid, fields, auxiliary_fields)
end

"""Evaluate one compiled Agate process diagnostic summed over selected tracers."""
@inline function process_tendency(
    plankton::NPDPlankton,
    diagnostics,
    tracers::Tuple,
    ::Val{Process},
    i, j, k, grid, fields, auxiliary_fields,
) where Process
    equations = getproperty(diagnostics, Process)
    tracer_values = _runtime_tracer_values(plankton, i, j, k, fields)
    auxiliary_values = _runtime_auxiliary_values(plankton, i, j, k, auxiliary_fields)
    t = zero(eltype(grid))
    return mapreduce(+, tracers; init=zero(t)) do tracer
        hasfield(typeof(equations), tracer) || return zero(t)
        x = zero(t)
        getfield(equations, tracer)(
            plankton.runtime, x, x, x, t, tracer_values..., auxiliary_values...
        )
    end
end
