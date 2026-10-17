using Adapt: adapt
import Adapt: adapt_structure

using ..ModelFamilies: AbstractModelFamily, plankton_roles
using ..Compilation: StaticZeroEquation
import ..Construction

import OceanBioME: chlorophyll
import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

using OceanBioME.Models.NutrientsPlanktonDetritusModels: NutrientsPlanktonDetritus
import OceanBioME.Models.NutrientsPlanktonDetritusModels:
    carbon_ratio,
    nitrogen_ratio,
    phosphate_ratio,
    iron_ratio,
    silicon_ratio,
    chlorophyll_ratio,
    dissolved_waste,
    inorganic_waste,
    nutrient_uptake,
    solid_waste
import OceanBioME.Models.NutrientsPlanktonDetritusModels.DetritusModels: grazing

"""OceanBioME NPD plankton component backed by one compiled Agate runtime.

`NPDPlankton` is an integration boundary, not a biological model. The compiled
Agate runtime owns the ecological equations; this wrapper maps their signed
tracer tendencies onto OceanBioME's existing NPD plankton hooks. Its `settings`
field holds the resolved Agate model settings used by the coupling.
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
    Settings,
}
    runtime::Runtime
    settings::Settings
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

function _validate_npd_settings(settings::NamedTuple)
    for name in (:carbon_ratio, :chlorophyll_ratio)
        hasproperty(settings, name) || throw(
            ArgumentError("NPDPlankton settings must define :$name."),
        )
        value = getproperty(settings, name)
        value isa Real && !(value isa Bool) && isfinite(value) && value >= 0 || throw(
            ArgumentError("NPDPlankton setting :$name must be finite and nonnegative."),
        )
    end
    settings.carbon_ratio > 0 || throw(ArgumentError("NPDPlankton carbon_ratio must be > 0."))
    return settings
end

# Resolve settings once on the host before constructing the statically routed adapter.
# Registered families carry resolved settings; standalone process definitions may not.
function _resolve_npd_settings(runtime, supplied::Union{Nothing,NamedTuple})
    resolved = runtime.metadata.model_settings
    if isempty(resolved)
        isnothing(supplied) && throw(ArgumentError(
            "NPDPlankton requires explicit settings when runtime.metadata.model_settings is empty.",
        ))
        return _validate_npd_settings(supplied)
    end

    _validate_npd_settings(resolved)
    isnothing(supplied) && return resolved
    _validate_npd_settings(supplied)
    names = keys(resolved)
    length(supplied) == length(resolved) &&
        all(name -> hasproperty(supplied, name), names) || throw(ArgumentError(
            "NPDPlankton settings keys must match runtime.metadata.model_settings keys $names.",
        ))
    for name in names
        actual = getproperty(supplied, name)
        expected = getproperty(resolved, name)
        # A Float64 keyword can round to the same Float32 scientific setting.
        equivalent = isequal(actual, expected) || (
            actual isa AbstractFloat && expected isa AbstractFloat &&
            isapprox(actual, expected;
                     rtol=max(eps(typeof(actual)), eps(typeof(expected))), atol=0)
        )
        equivalent || throw(ArgumentError(
            "NPDPlankton setting :$name must match runtime.metadata.model_settings (got $actual; expected $expected).",
        ))
    end
    return resolved
end

function _canonical_exchange_tracers(exchanges::NamedTuple)
    channels = (:solid, :dissolved, :inorganic)
    length(exchanges) == length(channels) &&
        all(channel -> hasproperty(exchanges, channel), channels) || throw(ArgumentError(
            "exchange_tracers must define exactly $channels; got $(keys(exchanges)).",
        ))
    return (
        solid=exchanges.solid,
        dissolved=exchanges.dissolved,
        inorganic=exchanges.inorganic,
    )
end

"""
    NPDPlankton(runtime; owned_components, phytoplankton_components=(), nutrient_tracers=(),
                exchange_tracers=_standard_exchange_tracers(runtime),
                consumed_detritus=(), dependencies=(), settings=nothing)

Wrap a compiled Agate runtime as an OceanBioME `NutrientsPlanktonDetritus` plankton component.
Component names are resolved once from Agate runtime metadata; all cell-level coupling is then
statically dispatched from the resulting tracer tuples. Every nonzero runtime tendency
must reach exactly one owned, nutrient, consumed-detritus or exchange hook. At NPD composition,
`nutrient_tracers` must match the native nutrient tracers and `consumed_detritus` must belong
to the native detritus component. Agate does not impose a nutrient-name vocabulary. A
dependency that also names a runtime auxiliary driver is read from the surrounding NPD tracer
fields rather than requested as a separate auxiliary field. Resolved family settings are
read from `runtime.metadata.model_settings`; any supplied `settings` must agree with them.
Standalone runtimes without resolved settings require explicit conversion ratios.
Exchange channel names are order-independent; an omitted channel is represented by `nothing`.
"""
function NPDPlankton(
    runtime;
    owned_components::Tuple,
    phytoplankton_components::Tuple=(),
    nutrient_tracers::Tuple=(),
    exchange_tracers::NamedTuple=_standard_exchange_tracers(runtime),
    consumed_detritus::Tuple=(),
    dependencies::Tuple=(),
    settings::Union{Nothing,NamedTuple}=nothing,
)
    exchange_tracers = _canonical_exchange_tracers(exchange_tracers)
    settings = _resolve_npd_settings(runtime, settings)
    owned = _component_tracers(runtime, owned_components)
    isempty(owned) && throw(ArgumentError("NPDPlankton must own at least one tracer."))
    phytoplankton = _component_tracers(runtime, phytoplankton_components)
    isempty(phytoplankton) && throw(
        ArgumentError("NPDPlankton must identify at least one phytoplankton tracer."),
    )
    runtime_tracers = required_biogeochemical_tracers(runtime)
    for (group, names) in (
        (:owned, owned), (:phytoplankton, phytoplankton),
        (:nutrient_tracers, nutrient_tracers), (:consumed_detritus, consumed_detritus),
        (:dependencies, dependencies),
    )
        all(name -> name isa Symbol, names) || throw(ArgumentError(
            "NPDPlankton $group must contain tracer names as Symbols.",
        ))
        length(unique(names)) == length(names) || throw(ArgumentError(
            "NPDPlankton $group contains duplicate names.",
        ))
    end
    all(name -> name in owned, phytoplankton) || throw(ArgumentError(
        "NPDPlankton phytoplankton tracers must belong to owned_components.",
    ))
    for tracer in nutrient_tracers
        tracer in runtime_tracers || throw(ArgumentError(
            "NPDPlankton nutrient tracer :$tracer is not an Agate runtime tracer.",
        ))
    end
    for tracer in consumed_detritus
        tracer in runtime_tracers || throw(ArgumentError(
            "NPDPlankton consumed detritus :$tracer is not an Agate runtime tracer.",
        ))
    end
    for (channel, tracer) in pairs(exchange_tracers)
        tracer === nothing && continue
        tracer isa Symbol && tracer in runtime_tracers || throw(ArgumentError(
            "NPDPlankton exchange channel :$channel names unknown runtime tracer :$tracer.",
        ))
    end
    exchanges = Tuple(values(exchange_tracers))
    # Each external tendency reaches OceanBioME through exactly one hook.
    routed = (
        owned..., nutrient_tracers..., consumed_detritus...,
        (name for name in exchanges if name !== nothing)...,
    )
    length(unique(routed)) == length(routed) || throw(ArgumentError(
        "NPDPlankton owned, nutrient, consumed-detritus and exchange tracers must not overlap.",
    ))
    # A nonzero Agate equation must have exactly one OceanBioME destination.
    # Static-zero runtime tracers may remain passive physical inputs.
    unrouted = Tuple(
        name for (name, equation) in pairs(runtime.equations)
        if !(equation isa StaticZeroEquation) && !(name in routed)
    )
    isempty(unrouted) || throw(ArgumentError(
        "NPDPlankton nonzero tendencies have no OceanBioME route: $unrouted.",
    ))
    owned_type = mapreduce(name -> typeof(Val(name)), (A, B) -> Union{A,B}, owned)

    return NPDPlankton{
        typeof(runtime),
        owned,
        owned_type,
        nutrient_tracers,
        exchanges,
        consumed_detritus,
        dependencies,
        phytoplankton,
        typeof(settings),
    }(runtime, settings)
end

"""Return family-specific `NPDPlankton` keyword overrides for a realized runtime.

Ownership, phytoplankton identity, external dependencies, standard exchange channels, and
resolved model settings are inferred from Agate metadata. External families only need to
declare choices that cannot be inferred safely, such as nutrient uptake tracers, consumed
detritus, or physical-driver dependencies. Configured `dependencies` are appended to the
inferred tracer dependencies.
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
    inferred = _npd_dependencies(runtime, configuration)
    configured = hasproperty(overrides, :dependencies) ? overrides.dependencies : ()
    configuration = merge(
        configuration, (; dependencies=_append_unique(inferred, configured))
    )
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

@inline _owned_tracers(::NPDPlankton{<:Any,OwnedTracers}) where {OwnedTracers} = OwnedTracers
@inline _nutrient_tracers(::NPDPlankton{<:Any,<:Any,<:Any,NutrientTracers}) where {NutrientTracers} = NutrientTracers
@inline _exchange_tracers(::NPDPlankton{<:Any,<:Any,<:Any,<:Any,ExchangeTracers}) where {ExchangeTracers} = ExchangeTracers
@inline _consumed_detritus(::NPDPlankton{<:Any,<:Any,<:Any,<:Any,<:Any,ConsumedDetritus}) where {ConsumedDetritus} = ConsumedDetritus
@inline _dependencies(::NPDPlankton{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,Dependencies}) where {Dependencies} = Dependencies
@inline phytoplankton_tracers(::NPDPlankton{<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,<:Any,PhytoplanktonTracers}) where {PhytoplanktonTracers} = PhytoplanktonTracers

@inline required_biogeochemical_tracers(plankton::NPDPlankton) = _owned_tracers(plankton)
@inline required_biogeochemical_auxiliary_fields(plankton::NPDPlankton) = Tuple(
    name for name in required_biogeochemical_auxiliary_fields(plankton.runtime)
    if !(name in _dependencies(plankton))
)
@inline biogeochemical_drift_velocity(plankton::NPDPlankton, tracer::Val) =
    biogeochemical_drift_velocity(plankton.runtime, tracer)

# OceanBioME's NPD drift fallback only checks detritus and inorganic carbon.
# Forward Agate-owned plankton tracers to the compiled runtime at the composed
# NPD boundary; other tracers retain OceanBioME's native drift dispatch.
@inline biogeochemical_drift_velocity(
    bgc::NutrientsPlanktonDetritus{<:Any,<:Any,PlanktonType},
    tracer::OwnedTracerType,
) where {
    OwnedTracerType<:Val,
    PlanktonType<:NPDPlankton{<:Any,<:Any,OwnedTracerType},
} = biogeochemical_drift_velocity(bgc.plankton, tracer)

@inline chlorophyll_ratio(plankton::NPDPlankton) = plankton.settings.chlorophyll_ratio
@inline carbon_ratio(
    plankton::NPDPlankton, ::NutrientsPlanktonDetritus{FloatType}
) where FloatType = convert(FloatType, plankton.settings.carbon_ratio)
@inline chlorophyll(plankton::NPDPlankton, model) = plankton.settings.chlorophyll_ratio *
    mapreduce(name -> getproperty(model.tracers, name), +, phytoplankton_tracers(plankton))

@inline function adapt_structure(
    to,
    plankton::NPDPlankton{
        <:Any,OwnedTracers,OwnedTracerType,NutrientTracers,ExchangeTracers,
        ConsumedDetritus,Dependencies,PhytoplanktonTracers,<:Any,
    },
) where {
    OwnedTracers,OwnedTracerType,NutrientTracers,ExchangeTracers,
    ConsumedDetritus,Dependencies,PhytoplanktonTracers,
}
    runtime = adapt(to, plankton.runtime)
    settings = adapt(to, plankton.settings)
    return NPDPlankton{
        typeof(runtime),OwnedTracers,OwnedTracerType,NutrientTracers,ExchangeTracers,
        ConsumedDetritus,Dependencies,PhytoplanktonTracers,typeof(settings),
    }(runtime, settings)
end

@inline function _append_unique(acc::Tuple, values::Tuple)
    isempty(values) && return acc
    head = first(values)
    next = head in acc ? acc : (acc..., head)
    return _append_unique(next, Base.tail(values))
end

# Validate the native component choices at the composed NPD boundary, where
# OceanBioME's nutrient and detritus tracer identities are finally known.
function _validate_npd_components(npd::NutrientsPlanktonDetritus, plankton::NPDPlankton)
    configured = _nutrient_tracers(plankton)
    nutrients = required_biogeochemical_tracers(npd.nutrients)
    length(configured) == length(nutrients) &&
        all(name -> name in nutrients, configured) &&
        all(name -> name in configured, nutrients) || throw(ArgumentError(
            "NPDPlankton nutrient_tracers $configured must match OceanBioME nutrient tracers $nutrients.",
        ))

    detritus = required_biogeochemical_tracers(npd.detritus)
    missing = Tuple(name for name in _consumed_detritus(plankton) if !(name in detritus))
    isempty(missing) || throw(ArgumentError(
        "NPDPlankton consumed_detritus $missing are absent from OceanBioME detritus tracers $detritus.",
    ))

    native = (
        nutrients..., detritus...,
        required_biogeochemical_tracers(npd.inorganic_carbon)...,
        required_biogeochemical_tracers(npd.oxygen)...,
    )
    collisions = Tuple(name for name in _owned_tracers(plankton) if name in native)
    isempty(collisions) || throw(ArgumentError(
        "NPDPlankton owned tracers $collisions also belong to native OceanBioME components.",
    ))
    return nothing
end

@inline function required_biogeochemical_tracers(
    npd::NutrientsPlanktonDetritus{<:Any,<:Any,PlanktonType},
) where {PlanktonType<:NPDPlankton}
    _validate_npd_components(npd, npd.plankton)
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

@inline function _runtime_auxiliary_value(
    ::Val{Auxiliary}, plankton::NPDPlankton, i, j, k, fields, auxiliary_fields
) where Auxiliary
    Auxiliary in _dependencies(plankton) &&
        return @inbounds getproperty(fields, Auxiliary)[i, j, k]
    return @inbounds getproperty(auxiliary_fields, Auxiliary)[i, j, k]
end

@inline function _runtime_auxiliary_values(
    plankton::NPDPlankton, i, j, k, fields, auxiliary_fields
)
    auxiliaries = required_biogeochemical_auxiliary_fields(plankton.runtime)
    return ntuple(Val(length(auxiliaries))) do n
        _runtime_auxiliary_value(
            Val(auxiliaries[n]), plankton, i, j, k, fields, auxiliary_fields
        )
    end
end

@inline function _agate_tendency(
    plankton::NPDPlankton, tracer::Val, i, j, k, t, fields, auxiliary_fields
)
    tracer_values = _runtime_tracer_values(plankton, i, j, k, fields)
    auxiliary_values = _runtime_auxiliary_values(
        plankton, i, j, k, fields, auxiliary_fields
    )
    x = zero(t)
    return plankton.runtime(tracer, x, x, x, t, tracer_values..., auxiliary_values...)
end

@inline _exchange_tendency(plankton, tracer, i, j, k, grid, fields, auxiliary_fields) =
    _agate_tendency(plankton, tracer, i, j, k, zero(eltype(grid)), fields, auxiliary_fields)

# Restrict this overload to Agate-owned living tracers. The Val bound makes it
# more specific than OceanBioME's generic Val fallback, while native nutrient,
# detritus, carbon, and oxygen tracer dispatch remains intact.
@inline (bgc::NutrientsPlanktonDetritus{<:Any,<:Any,PlanktonType})(
    i, j, k, grid, tracer::OwnedTracerType, clock, fields, auxiliary_fields
) where {
    OwnedTracerType<:Val,
    PlanktonType<:NPDPlankton{<:Any,<:Any,OwnedTracerType},
} = _agate_tendency(
    bgc.plankton, tracer, i, j, k, clock.time, fields, auxiliary_fields
)

@inline function _npd_nutrient_uptake(
    i, j, k, grid, nutrient::Val{Nutrient}, plankton::NPDPlankton,
    fields, auxiliary_fields,
) where Nutrient
    Nutrient in _nutrient_tracers(plankton) || return zero(eltype(grid))
    return -_exchange_tendency(plankton, nutrient, i, j, k, grid, fields, auxiliary_fields)
end

@inline nutrient_uptake(
    i, j, k, grid, nutrient::Val{Nutrient}, plankton::NPDPlankton,
    ::NutrientsPlanktonDetritus, fields, auxiliary_fields,
) where Nutrient = _npd_nutrient_uptake(
    i, j, k, grid, nutrient, plankton, fields, auxiliary_fields
)

# OceanBioME has Val-specific methods for elemental nutrients. Resolve their
# intersection with NPDPlankton, preserving OceanBioME's stoichiometric fallback
# when an element is not itself an Agate nutrient tracer (notably carbon).
for (nutrient, ratio) in (
    (:N, :nitrogen_ratio), (:PO₄, :phosphate_ratio), (:Fe, :iron_ratio),
    (:Si, :silicon_ratio), (:C, :carbon_ratio),
)
    @eval @inline function nutrient_uptake(
        i, j, k, grid, ::Val{$(QuoteNode(nutrient))}, plankton::NPDPlankton,
        bgc::NutrientsPlanktonDetritus, fields, auxiliary_fields,
    )
        $(QuoteNode(nutrient)) in _nutrient_tracers(plankton) &&
            return _npd_nutrient_uptake(
                i, j, k, grid, Val($(QuoteNode(nutrient))), plankton, fields, auxiliary_fields
            )
        return $ratio(i, j, k, grid, plankton, bgc, fields) *
            nutrient_uptake(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
    end
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
