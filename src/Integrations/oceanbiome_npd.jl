using Adapt: adapt
import Adapt: adapt_structure

using ..ModelFamilies: AbstractModelFamily, plankton_roles
import ..Construction

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
                consumed_detritus=(), dependencies=(), traits)

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
    runtime_tracers = required_biogeochemical_tracers(runtime)
    for (channel, tracer) in pairs(exchange_tracers)
        tracer === nothing && continue
        tracer in runtime_tracers || throw(ArgumentError(
            "NPDPlankton exchange channel :$channel names unknown runtime tracer :$tracer.",
        ))
    end
    exchanges = Tuple(values(exchange_tracers))
    _validate_npd_traits(traits)

    return NPDPlankton{
        typeof(runtime),
        owned,
        owned_type,
        nutrient_tracers,
        exchanges,
        consumed_detritus,
        dependencies,
        phytoplankton,
        typeof(traits),
    }(runtime, traits)
end

"""Return family-specific `NPDPlankton` keyword overrides for a realized runtime.

Ownership, phytoplankton identity, external dependencies, standard exchange channels, and
resolved model settings are inferred from Agate metadata. External families only need to
declare choices that cannot be inferred safely, such as nutrient uptake tracers or
consumed detritus.
"""
npd_configuration(::AbstractModelFamily, _runtime) = (;)

function _standard_exchange_tracers(runtime)
    tracers = required_biogeochemical_tracers(runtime)
    present(name) = name in tracers ? name : nothing
    return (
        solid=present(:solid_waste),
        dissolved=present(:dissolved_waste),
        inorganic=present(:inorganic_waste),
    )
end

function _npd_default_configuration(family::AbstractModelFamily, runtime)
    roles = plankton_roles(family)
    hasproperty(roles, :phytoplankton) || throw(
        ArgumentError("NPD plankton families must define a :phytoplankton role."),
    )
    return (;
        owned_components=Tuple(unique(values(roles))),
        phytoplankton_components=(roles.phytoplankton,),
        nutrient_tracers=(),
        exchange_tracers=_standard_exchange_tracers(runtime),
        consumed_detritus=(),
        traits=runtime.metadata.model_settings,
    )
end

function _npd_dependencies(runtime, configuration)
    owned = _component_tracers(runtime, configuration.owned_components)
    exchanges = Tuple(
        value for value in values(configuration.exchange_tracers) if value !== nothing
    )
    return Tuple(
        tracer for tracer in required_biogeochemical_tracers(runtime)
        if !(tracer in owned) && !(tracer in exchanges)
    )
end

function _wrap_npd_runtime(family::AbstractModelFamily, runtime)
    overrides = npd_configuration(family, runtime)
    configuration = merge(_npd_default_configuration(family, runtime), overrides)
    if !hasproperty(overrides, :dependencies)
        configuration = merge(
            configuration, (; dependencies=_npd_dependencies(runtime, configuration))
        )
    end
    return NPDPlankton(runtime; configuration...)
end

_require_npd_grid(sinking_tracers, grid) =
    !isnothing(sinking_tracers) && isnothing(grid) ?
        throw(ArgumentError("grid is required when `sinking_tracers` are configured")) : nothing

"""Construct a registered Agate family directly as an OceanBioME NPD plankton component."""
function construct_npd_plankton(
    family::AbstractModelFamily;
    plankton_pfts::NamedTuple,
    parameter_overrides::NamedTuple=(;),
    setting_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
)
    _require_npd_grid(sinking_tracers, grid)
    runtime = Construction.construct(
        family;
        plankton_pfts,
        parameter_overrides,
        setting_overrides,
        sinking_tracers,
        open_bottom,
        grid,
        arch,
        scalar_type,
    )
    return _wrap_npd_runtime(family, runtime)
end

"""Construct an NPD plankton component and capture the canonical Agate family recipe."""
function construct_npd_plankton_plus_recipe(
    family::AbstractModelFamily;
    plankton_pfts::NamedTuple,
    parameter_overrides::NamedTuple=(;),
    setting_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    grid=nothing,
    arch=nothing,
    scalar_type=nothing,
)
    _require_npd_grid(sinking_tracers, grid)
    runtime, recipe = Construction.construct_plus_recipe(
        family;
        plankton_pfts,
        parameter_overrides,
        setting_overrides,
        sinking_tracers,
        open_bottom,
        grid,
        arch,
        scalar_type,
    )
    return _wrap_npd_runtime(family, runtime), recipe
end

"""Replay a registered Agate family recipe directly as an OceanBioME NPD plankton component."""
function construct_npd_plankton(
    recipe::Construction.ModelRecipe; grid=nothing, arch=nothing, scalar_type=nothing
)
    family = Construction.replay_family(recipe)
    _require_npd_grid(recipe.sinking_tracers, grid)
    runtime = Construction.construct(
        recipe;
        grid,
        arch,
        scalar_type,
    )
    return _wrap_npd_runtime(family, runtime)
end

@inline _owned_tracers(::NPDPlankton{R,O}) where {R,O} = O
@inline _nutrient_tracers(::NPDPlankton{R,O,OT,N}) where {R,O,OT,N} = N
@inline _exchange_tracers(::NPDPlankton{R,O,OT,N,E}) where {R,O,OT,N,E} = E
@inline _consumed_detritus(::NPDPlankton{R,O,OT,N,E,D}) where {R,O,OT,N,E,D} = D
@inline _dependencies(::NPDPlankton{R,O,OT,N,E,D,Deps}) where {R,O,OT,N,E,D,Deps} = Deps
@inline phytoplankton_tracers(::NPDPlankton{R,O,OT,N,E,D,Deps,P}) where {R,O,OT,N,E,D,Deps,P} = P

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

@inline function adapt_structure(to, plankton::NPDPlankton{R,O,OT,N,E,D,Deps,P,T}) where {R,O,OT,N,E,D,Deps,P,T}
    runtime = adapt(to, plankton.runtime)
    traits = adapt(to, plankton.traits)
    return NPDPlankton{typeof(runtime),O,OT,N,E,D,Deps,P,typeof(traits)}(runtime, traits)
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
    Runtime,OwnedTracers,OwnedTracerType,N,E,D,Deps,P,T,
    PLA<:NPDPlankton{Runtime,OwnedTracers,OwnedTracerType,N,E,D,Deps,P,T},
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
