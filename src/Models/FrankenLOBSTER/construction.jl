using ...Construction
import ...Integrations

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

Integrations.npd_diagnostic_processes(::FrankenLOBSTERFamily) = _CALCITE_DIAGNOSTIC_PROCESSES

function Integrations.npd_configuration(
    ::FrankenLOBSTERFamily, runtime, parameters::NamedTuple
)
    return (;
        owned_components=(:P, :Z, :H),
        phytoplankton_components=(:P,),
        nutrient_tracers=(:NO₃, :NH₄),
        exchange_tracers=(
            solid=:solid_waste, dissolved=:dissolved_waste, inorganic=:inorganic_waste,
        ),
        consumed_detritus=(:DOM,),
        dependencies=(:NO₃, :NH₄, :DOM, :T),
        traits=_traits(parameters),
        coupling=FrankenLOBSTERCoupling(runtime.metadata.process_diagnostics),
    )
end

function _construction_inputs(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    parameters::NamedTuple=(;),
    grid=nothing,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    family = FrankenLOBSTERFamily()
    realization = (;
        plankton_pfts=Construction.plankton_realization(family, size_structure),
        parameter_overrides=parameters,
        sinking_tracers,
        open_bottom,
    )
    return (; family, realization, execution=(; grid))
end

"""Construct the Agate plankton component for composition with OceanBioME `LOBSTER`."""
function construct(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    return Integrations.construct_npd_plankton(
        inputs.family; inputs.realization..., inputs.execution...
    )
end

"""Construct FrankenLOBSTER plankton and capture its versioned Agate recipe."""
function construct_plus_recipe(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    return Integrations.construct_npd_plankton_plus_recipe(
        inputs.family; inputs.realization..., inputs.execution...
    )
end

"""Replay a FrankenLOBSTER recipe into the Agate plankton component."""
function construct(recipe::Construction.ModelRecipe; grid=nothing)
    recipe.family == :FrankenLOBSTER || throw(ArgumentError(
        "FrankenLOBSTER.construct requires a FrankenLOBSTER recipe; got $(recipe.family)"
    ))
    return Integrations.construct_npd_plankton(recipe; grid)
end
