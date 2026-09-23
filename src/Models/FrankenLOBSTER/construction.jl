using OceanBioME: BoxModelGrid
using OceanBioME.Models.NutrientsPlanktonDetritusModels: LOBSTER
using ...Construction

const _SIZE_ROLES = (:phytoplankton, :zooplankton, :bacterioplankton)

const _CALCITE_DIAGNOSTIC_PROCESSES = (
    :nitrate_growth_P, :ammonium_growth_P, :grazing_Z_on_living, :mortality_P,
)

function _validate_nonnegative_finite(name, value)
    value isa Real && !(value isa Bool) && isfinite(value) && value >= zero(value) ||
        throw(ArgumentError("$name must be a finite nonnegative real number; got $(repr(value))"))
    return value
end

function _validate_fraction(name, value)
    _validate_nonnegative_finite(name, value)
    value <= one(value) || throw(ArgumentError("$name must be <= 1; got $(repr(value))"))
    return value
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

    # Generic Agate realization validates non-empty roles, duplicate PFT identities, and
    # SizeClass specifications after this user-facing role -> logical component mapping.
    return (P=normalize(:phytoplankton), Z=normalize(:zooplankton), H=normalize(:bacterioplankton))
end

function _construct_plankton(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    grid=BoxModelGrid(),
    parameters::NamedTuple=(;),
    phytoplankton_chlorophyll_ratio=1.31,
    carbon_ratio=FRANKENLOBSTER_CARBON_RATIO,
    calcium_carbonate_rain_ratio=FRANKENLOBSTER_CALCIUM_CARBONATE_RAIN_RATIO,
    zooplankton_calcium_carbonate_dissolution=FRANKENLOBSTER_ZOOPLANKTON_CALCIUM_CARBONATE_DISSOLUTION,
    scalar_type=nothing,
    arch=nothing,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    realization = _plankton_realization(size_structure)
    _validate_nonnegative_finite("carbon_ratio", carbon_ratio)
    carbon_ratio > zero(carbon_ratio) || throw(ArgumentError("carbon_ratio must be > 0"))
    _validate_nonnegative_finite("calcium_carbonate_rain_ratio", calcium_carbonate_rain_ratio)
    _validate_fraction(
        "zooplankton_calcium_carbonate_dissolution", zooplankton_calcium_carbonate_dissolution
    )
    runtime = Construction.construct(
        FrankenLOBSTERFamily();
        plankton_pfts=realization,
        grid,
        parameter_overrides=parameters,
        sinking_tracers,
        open_bottom,
        scalar_type,
        arch,
        diagnostic_processes=_CALCITE_DIAGNOSTIC_PROCESSES,
    )
    phytoplankton_tracers = Tuple(
        tracer
        for pft in keys(realization.P)
        for tracer in getproperty(runtime.metadata.pft_entities, pft)
    )
    return FrankenLOBSTERPlankton(
        runtime,
        runtime.metadata.plankton_tracers,
        (:solid_waste, :inorganic_waste, :dissolved_waste);
        phytoplankton_tracers,
        process_diagnostics=runtime.metadata.process_diagnostics,
        chlorophyll_ratio=convert(eltype(grid), phytoplankton_chlorophyll_ratio),
        carbon_ratio=convert(eltype(grid), carbon_ratio),
        calcium_carbonate_rain_ratio=convert(eltype(grid), calcium_carbonate_rain_ratio),
        zooplankton_calcium_carbonate_dissolution=convert(
            eltype(grid), zooplankton_calcium_carbonate_dissolution
        ),
    )
end

"""Construct coupled FrankenLOBSTER: Agate living ecology inside OceanBioME LOBSTER."""
function construct(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    parameters::NamedTuple=(;),
    phytoplankton_chlorophyll_ratio=1.31,
    carbon_ratio=FRANKENLOBSTER_CARBON_RATIO,
    calcium_carbonate_rain_ratio=FRANKENLOBSTER_CALCIUM_CARBONATE_RAIN_RATIO,
    zooplankton_calcium_carbonate_dissolution=FRANKENLOBSTER_ZOOPLANKTON_CALCIUM_CARBONATE_DISSOLUTION,
    grid=BoxModelGrid(),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    kwargs...,
)
    plankton = _construct_plankton(;
        size_structure, parameters, phytoplankton_chlorophyll_ratio, carbon_ratio,
        calcium_carbonate_rain_ratio, zooplankton_calcium_carbonate_dissolution, grid,
        sinking_tracers, open_bottom,
    )
    return LOBSTER(
        grid;
        limiting_nutrients=(:nitrate, :ammonia, :iron),
        plankton,
        open_bottom,
        kwargs...,
    )
end
