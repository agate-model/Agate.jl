using ...Construction
using ...Integrations: NPDPlankton

const _CALCITE_DIAGNOSTIC_PROCESSES = (
    :nitrate_growth_P, :ammonia_growth_P, :grazing_Z_on_living, :mortality_P,
)
const _TRAIT_DEFAULTS = (
    chlorophyll_ratio=FRANKENLOBSTER_CHLOROPHYLL_RATIO,
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

_require_sinking_grid(sinking_tracers, grid) =
    !isnothing(sinking_tracers) && isnothing(grid) ?
        throw(ArgumentError("grid is required when `sinking_tracers` are configured")) : nothing

function _wrap_plankton(runtime, parameters, ::Type{T}) where T
    traits = _traits(parameters)
    typed_traits = NamedTuple{keys(traits)}(map(value -> convert(T, value), values(traits)))
    return NPDPlankton(
        runtime;
        owned_components=(:P, :Z, :H),
        phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃, :NH₄),
        exchange_tracers=(
            solid=:solid_waste, dissolved=:dissolved_waste, inorganic=:inorganic_waste,
        ),
        consumed_detritus=(:DOM,),
        dependencies=(:NO₃, :NH₄, :DOM, :T),
        traits=typed_traits,
        coupling=FrankenLOBSTERCoupling(runtime.metadata.process_diagnostics),
    )
end

function _construct_plankton(realization, parameters; grid=nothing, sinking_tracers=nothing, open_bottom=true)
    _require_sinking_grid(sinking_tracers, grid)
    runtime = Construction.construct(
        FrankenLOBSTERFamily(); plankton_pfts=realization, grid,
        parameter_overrides=_without_traits(parameters), sinking_tracers, open_bottom,
        diagnostic_processes=_CALCITE_DIAGNOSTIC_PROCESSES,
    )
    return _wrap_plankton(runtime, parameters, isnothing(grid) ? Float64 : eltype(grid))
end

"""Construct the Agate plankton component for composition with OceanBioME `LOBSTER`."""
function construct(; size_structure=DEFAULT_SIZE_STRUCTURE, parameters::NamedTuple=(;),
                   grid=nothing, sinking_tracers=nothing, open_bottom::Bool=true)
    return _construct_plankton(
        Construction.plankton_realization(FrankenLOBSTERFamily(), size_structure), parameters;
        grid, sinking_tracers, open_bottom
    )
end

"""Construct FrankenLOBSTER plankton and capture its versioned Agate recipe."""
function construct_plus_recipe(; size_structure=DEFAULT_SIZE_STRUCTURE, parameters::NamedTuple=(;),
                               grid=nothing, sinking_tracers=nothing, open_bottom::Bool=true)
    realization = Construction.plankton_realization(FrankenLOBSTERFamily(), size_structure)
    _require_sinking_grid(sinking_tracers, grid)
    runtime, recipe = Construction.construct_plus_recipe(
        FrankenLOBSTERFamily(); plankton_pfts=realization, parameter_overrides=parameters,
        sinking_tracers, open_bottom, grid, diagnostic_processes=_CALCITE_DIAGNOSTIC_PROCESSES,
    )
    return _wrap_plankton(
        runtime, parameters, isnothing(grid) ? Float64 : eltype(grid)
    ), recipe
end

"""Replay a FrankenLOBSTER recipe into the Agate plankton component."""
function construct(recipe::Construction.ModelRecipe; grid=nothing)
    recipe.family == :FrankenLOBSTER || throw(ArgumentError(
        "FrankenLOBSTER.construct requires a FrankenLOBSTER recipe; got $(recipe.family)"
    ))
    _require_sinking_grid(recipe.sinking_tracers, grid)
    runtime = Construction.construct(recipe; grid, diagnostic_processes=_CALCITE_DIAGNOSTIC_PROCESSES)
    return _wrap_plankton(
        runtime, recipe.parameter_overrides, isnothing(grid) ? Float64 : eltype(grid)
    )
end
