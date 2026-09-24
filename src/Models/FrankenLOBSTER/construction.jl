using ...Construction

const _SIZE_ROLES = (:phytoplankton, :zooplankton, :bacterioplankton)
const _CALCITE_DIAGNOSTIC_PROCESSES = (
    :nitrate_growth_P, :ammonia_growth_P, :grazing_Z_on_living, :mortality_P,
)
const _TRAIT_DEFAULTS = (
    phytoplankton_chlorophyll_ratio=FRANKENLOBSTER_CHLOROPHYLL_RATIO,
    carbon_ratio=FRANKENLOBSTER_CARBON_RATIO,
    calcium_carbonate_rain_ratio=FRANKENLOBSTER_CALCIUM_CARBONATE_RAIN_RATIO,
    zooplankton_calcium_carbonate_dissolution=FRANKENLOBSTER_ZOOPLANKTON_CALCIUM_CARBONATE_DISSOLUTION,
)
const _TRAIT_NAMES = keys(_TRAIT_DEFAULTS)

function _without_traits(parameters::NamedTuple)
    names = Tuple(name for name in keys(parameters) if !(name in _TRAIT_NAMES))
    return NamedTuple{names}(Tuple(getproperty(parameters, name) for name in names))
end
Construction.recipe_runtime_parameter_overrides(::FrankenLOBSTERFamily, overrides::NamedTuple) =
    _without_traits(overrides)

function _traits(parameters::NamedTuple)
    merged = merge(_TRAIT_DEFAULTS, parameters)
    traits = NamedTuple{_TRAIT_NAMES}(Tuple(getproperty(merged, name) for name in _TRAIT_NAMES))
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && x >= 0, values(traits)) ||
        throw(ArgumentError("FrankenLOBSTER traits must be finite and nonnegative"))
    traits.carbon_ratio > 0 || throw(ArgumentError("carbon_ratio must be > 0"))
    traits.zooplankton_calcium_carbonate_dissolution <= 1 ||
        throw(ArgumentError("zooplankton_calcium_carbonate_dissolution must be <= 1"))
    return traits
end

function _plankton_realization(size_structure)
    size_structure isa NamedTuple || throw(ArgumentError("size_structure must be a NamedTuple"))
    Set(keys(size_structure)) == Set(_SIZE_ROLES) || throw(ArgumentError(
        "size_structure must define exactly phytoplankton, zooplankton, and bacterioplankton"
    ))
    normalize(role) = begin
        pfts = getproperty(size_structure, role)
        pfts isa NamedTuple || throw(ArgumentError("size_structure.$role must be a NamedTuple"))
        NamedTuple{keys(pfts)}(map(Construction.normalize_pft_size_structure, values(pfts)))
    end
    return (P=normalize(:phytoplankton), Z=normalize(:zooplankton), H=normalize(:bacterioplankton))
end

_require_sinking_grid(sinking_tracers, grid) =
    !isnothing(sinking_tracers) && isnothing(grid) ?
        throw(ArgumentError("grid is required when `sinking_tracers` are configured")) : nothing

function _wrap_plankton(runtime, realization, parameters, ::Type{T}) where T
    traits = _traits(parameters)
    typed_traits = NamedTuple{keys(traits)}(map(value -> convert(T, value), values(traits)))
    phytoplankton_tracers = Tuple(
        tracer for pft in keys(realization.P) for tracer in getproperty(runtime.metadata.pft_entities, pft)
    )
    return FrankenLOBSTERPlankton(
        runtime, runtime.metadata.plankton_tracers, (:solid_waste, :inorganic_waste, :dissolved_waste);
        phytoplankton_tracers, process_diagnostics=runtime.metadata.process_diagnostics, traits=typed_traits,
    )
end

function _construct_plankton(realization, parameters; grid=nothing, sinking_tracers=nothing, open_bottom=true)
    _require_sinking_grid(sinking_tracers, grid)
    runtime = Construction.construct(
        FrankenLOBSTERFamily(); plankton_pfts=realization, grid,
        parameter_overrides=_without_traits(parameters), sinking_tracers, open_bottom,
        diagnostic_processes=_CALCITE_DIAGNOSTIC_PROCESSES,
    )
    return _wrap_plankton(runtime, realization, parameters, isnothing(grid) ? Float64 : eltype(grid))
end

"""Construct the Agate plankton component for composition with OceanBioME `LOBSTER`."""
function construct(; size_structure=DEFAULT_SIZE_STRUCTURE, parameters::NamedTuple=(;),
                   grid=nothing, sinking_tracers=nothing, open_bottom::Bool=true)
    return _construct_plankton(
        _plankton_realization(size_structure), parameters; grid, sinking_tracers, open_bottom
    )
end

"""Construct FrankenLOBSTER plankton and capture its versioned Agate recipe."""
function construct_plus_recipe(; size_structure=DEFAULT_SIZE_STRUCTURE, parameters::NamedTuple=(;),
                               grid=nothing, sinking_tracers=nothing, open_bottom::Bool=true)
    realization = _plankton_realization(size_structure)
    recipe = Construction.capture_model_recipe(
        FrankenLOBSTERFamily(); plankton_pfts=realization,
        parameter_overrides=parameters, sinking_tracers, open_bottom,
    )
    return _construct_plankton(realization, parameters; grid, sinking_tracers, open_bottom), recipe
end

"""Replay a FrankenLOBSTER recipe into the Agate plankton component."""
function construct(recipe::Construction.ModelRecipe; grid=nothing)
    recipe.family == :FrankenLOBSTER || throw(ArgumentError(
        "FrankenLOBSTER.construct requires a FrankenLOBSTER recipe; got $(recipe.family)"
    ))
    _require_sinking_grid(recipe.sinking_tracers, grid)
    runtime = Construction.construct(recipe; grid, diagnostic_processes=_CALCITE_DIAGNOSTIC_PROCESSES)
    return _wrap_plankton(
        runtime, recipe.plankton_pfts, recipe.parameter_overrides, isnothing(grid) ? Float64 : eltype(grid)
    )
end
