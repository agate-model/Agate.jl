using OceanBioME: BoxModelGrid
using OceanBioME.Models.NutrientsPlanktonDetritusModels: LOBSTER
using ...Construction

const _SIZE_ROLES = (:phytoplankton, :zooplankton, :bacterioplankton)
const _CALCITE_DIAGNOSTIC_PROCESSES = (
    :nitrate_growth_P, :ammonium_growth_P, :grazing_Z_on_living, :mortality_P,
)
const _COUPLING_PARAMETER_DEFAULTS = (
    phytoplankton_chlorophyll_ratio=FRANKENLOBSTER_CHLOROPHYLL_RATIO,
    carbon_ratio=FRANKENLOBSTER_CARBON_RATIO,
    calcium_carbonate_rain_ratio=FRANKENLOBSTER_CALCIUM_CARBONATE_RAIN_RATIO,
    zooplankton_calcium_carbonate_dissolution=
        FRANKENLOBSTER_ZOOPLANKTON_CALCIUM_CARBONATE_DISSOLUTION,
)

function _validate_nonnegative_finite(name, value)
    value isa Real && !(value isa Bool) && isfinite(value) && value >= zero(value) ||
        throw(ArgumentError("$name must be a finite nonnegative real number; got $(repr(value))"))
    return value
end

function _runtime_parameter_overrides(parameters::NamedTuple)
    coupling_names = keys(_COUPLING_PARAMETER_DEFAULTS)
    names = Tuple(name for name in keys(parameters) if !(name in coupling_names))
    return NamedTuple{names}(Tuple(getproperty(parameters, name) for name in names))
end

Construction.recipe_runtime_parameter_overrides(
    ::FrankenLOBSTERFamily, overrides::NamedTuple
) = _runtime_parameter_overrides(overrides)

function _coupling_values(parameters::NamedTuple)
    values = merge(_COUPLING_PARAMETER_DEFAULTS, parameters)
    names = keys(_COUPLING_PARAMETER_DEFAULTS)
    coupling = NamedTuple{names}(Tuple(getproperty(values, name) for name in names))
    _validate_nonnegative_finite(
        "phytoplankton_chlorophyll_ratio", coupling.phytoplankton_chlorophyll_ratio
    )
    _validate_nonnegative_finite("carbon_ratio", coupling.carbon_ratio)
    coupling.carbon_ratio > zero(coupling.carbon_ratio) ||
        throw(ArgumentError("carbon_ratio must be > 0"))
    _validate_nonnegative_finite(
        "calcium_carbonate_rain_ratio", coupling.calcium_carbonate_rain_ratio
    )
    fraction = coupling.zooplankton_calcium_carbonate_dissolution
    _validate_nonnegative_finite("zooplankton_calcium_carbonate_dissolution", fraction)
    fraction <= one(fraction) || throw(ArgumentError(
        "zooplankton_calcium_carbonate_dissolution must be <= 1; got $(repr(fraction))"
    ))
    return coupling
end

function _parameter_overrides(
    parameters::NamedTuple;
    phytoplankton_chlorophyll_ratio=FRANKENLOBSTER_CHLOROPHYLL_RATIO,
    carbon_ratio=FRANKENLOBSTER_CARBON_RATIO,
    calcium_carbonate_rain_ratio=FRANKENLOBSTER_CALCIUM_CARBONATE_RAIN_RATIO,
    zooplankton_calcium_carbonate_dissolution=
        FRANKENLOBSTER_ZOOPLANKTON_CALCIUM_CARBONATE_DISSOLUTION,
)
    keyword_values = (;
        phytoplankton_chlorophyll_ratio, carbon_ratio, calcium_carbonate_rain_ratio,
        zooplankton_calcium_carbonate_dissolution,
    )
    overrides = parameters
    for name in keys(keyword_values)
        value = getproperty(keyword_values, name)
        value == getproperty(_COUPLING_PARAMETER_DEFAULTS, name) && continue
        hasproperty(overrides, name) && throw(ArgumentError(
            "parameter :$name cannot be supplied through both `parameters` and `$name`"
        ))
        overrides = merge(overrides, NamedTuple{(name,)}((value,)))
    end
    _coupling_values(overrides)
    return overrides
end

function _plankton_realization(size_structure)
    size_structure isa NamedTuple || throw(ArgumentError("size_structure must be a NamedTuple"))
    Set(keys(size_structure)) == Set(_SIZE_ROLES) || throw(ArgumentError(
        "size_structure must define exactly phytoplankton, zooplankton, and bacterioplankton"
    ))

    normalize(role) = begin
        pfts = getproperty(size_structure, role)
        pfts isa NamedTuple || throw(ArgumentError("size_structure.$role must be a NamedTuple"))
        NamedTuple{keys(pfts)}(Tuple(
            Construction.normalize_pft_size_structure(value) for value in values(pfts)
        ))
    end
    return (P=normalize(:phytoplankton), Z=normalize(:zooplankton), H=normalize(:bacterioplankton))
end

function _wrap_plankton(runtime, realization, parameters, grid)
    coupling = _coupling_values(parameters)
    phytoplankton_tracers = Tuple(
        tracer for pft in keys(realization.P)
        for tracer in getproperty(runtime.metadata.pft_entities, pft)
    )
    return FrankenLOBSTERPlankton(
        runtime,
        runtime.metadata.plankton_tracers,
        (:solid_waste, :inorganic_waste, :dissolved_waste);
        phytoplankton_tracers,
        process_diagnostics=runtime.metadata.process_diagnostics,
        chlorophyll_ratio=convert(eltype(grid), coupling.phytoplankton_chlorophyll_ratio),
        carbon_ratio=convert(eltype(grid), coupling.carbon_ratio),
        calcium_carbonate_rain_ratio=convert(eltype(grid), coupling.calcium_carbonate_rain_ratio),
        zooplankton_calcium_carbonate_dissolution=convert(
            eltype(grid), coupling.zooplankton_calcium_carbonate_dissolution
        ),
    )
end

function _construct_plankton(
    realization, parameters, grid; sinking_tracers=nothing, open_bottom=true
)
    runtime = Construction.construct(
        FrankenLOBSTERFamily();
        plankton_pfts=realization,
        grid,
        parameter_overrides=_runtime_parameter_overrides(parameters),
        sinking_tracers,
        open_bottom,
        diagnostic_processes=_CALCITE_DIAGNOSTIC_PROCESSES,
    )
    return _wrap_plankton(runtime, realization, parameters, grid)
end

function _coupled(plankton, grid, open_bottom; kwargs...)
    return LOBSTER(
        grid;
        limiting_nutrients=(:nitrate, :ammonia, :iron),
        plankton,
        open_bottom,
        kwargs...,
    )
end

function _inputs(
    size_structure, parameters, sinking_tracers, open_bottom;
    phytoplankton_chlorophyll_ratio, carbon_ratio, calcium_carbonate_rain_ratio,
    zooplankton_calcium_carbonate_dissolution,
)
    parameters = _parameter_overrides(
        parameters;
        phytoplankton_chlorophyll_ratio, carbon_ratio, calcium_carbonate_rain_ratio,
        zooplankton_calcium_carbonate_dissolution,
    )
    realization = _plankton_realization(size_structure)
    return realization, parameters, (; plankton_pfts=realization, parameter_overrides=parameters,
                                      sinking_tracers, open_bottom)
end

"""Construct coupled FrankenLOBSTER: Agate living ecology inside OceanBioME LOBSTER."""
function construct(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    parameters::NamedTuple=(;),
    phytoplankton_chlorophyll_ratio=FRANKENLOBSTER_CHLOROPHYLL_RATIO,
    carbon_ratio=FRANKENLOBSTER_CARBON_RATIO,
    calcium_carbonate_rain_ratio=FRANKENLOBSTER_CALCIUM_CARBONATE_RAIN_RATIO,
    zooplankton_calcium_carbonate_dissolution=
        FRANKENLOBSTER_ZOOPLANKTON_CALCIUM_CARBONATE_DISSOLUTION,
    grid=BoxModelGrid(),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    kwargs...,
)
    realization, parameters, _ = _inputs(
        size_structure, parameters, sinking_tracers, open_bottom;
        phytoplankton_chlorophyll_ratio, carbon_ratio, calcium_carbonate_rain_ratio,
        zooplankton_calcium_carbonate_dissolution,
    )
    plankton = _construct_plankton(realization, parameters, grid; sinking_tracers, open_bottom)
    return _coupled(plankton, grid, open_bottom; kwargs...)
end

"""Construct FrankenLOBSTER and return the coupled model with its versioned family recipe."""
function construct_plus_recipe(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    parameters::NamedTuple=(;),
    phytoplankton_chlorophyll_ratio=FRANKENLOBSTER_CHLOROPHYLL_RATIO,
    carbon_ratio=FRANKENLOBSTER_CARBON_RATIO,
    calcium_carbonate_rain_ratio=FRANKENLOBSTER_CALCIUM_CARBONATE_RAIN_RATIO,
    zooplankton_calcium_carbonate_dissolution=
        FRANKENLOBSTER_ZOOPLANKTON_CALCIUM_CARBONATE_DISSOLUTION,
    grid=BoxModelGrid(),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    kwargs...,
)
    realization, parameters, recipe_inputs = _inputs(
        size_structure, parameters, sinking_tracers, open_bottom;
        phytoplankton_chlorophyll_ratio, carbon_ratio, calcium_carbonate_rain_ratio,
        zooplankton_calcium_carbonate_dissolution,
    )
    recipe = Construction.capture_model_recipe(FrankenLOBSTERFamily(); recipe_inputs...)
    plankton = _construct_plankton(realization, parameters, grid; sinking_tracers, open_bottom)
    return _coupled(plankton, grid, open_bottom; kwargs...), recipe
end

"""Replay a FrankenLOBSTER recipe with caller-owned OceanBioME coupling configuration."""
function construct_from_recipe(recipe::Construction.ModelRecipe; grid=BoxModelGrid(), kwargs...)
    recipe.family == :FrankenLOBSTER || throw(ArgumentError(
        "FrankenLOBSTER.construct_from_recipe requires a FrankenLOBSTER recipe; " *
        "got $(recipe.family)"
    ))
    runtime = Construction.construct(
        recipe; grid, diagnostic_processes=_CALCITE_DIAGNOSTIC_PROCESSES
    )
    plankton = _wrap_plankton(runtime, recipe.plankton_pfts, recipe.parameter_overrides, grid)
    return _coupled(plankton, grid, recipe.open_bottom; kwargs...)
end
