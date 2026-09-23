using OceanBioME: BoxModelGrid
using OceanBioME.Models.NutrientsPlanktonDetritusModels: LOBSTER
using ...Construction

const _SIZE_ROLES = (:phytoplankton, :zooplankton, :bacterioplankton)

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
    scalar_type=nothing,
    arch=nothing,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    realization = _plankton_realization(size_structure)
    runtime = Construction.construct(
        FrankenLOBSTERFamily();
        plankton_pfts=realization,
        grid,
        parameter_overrides=parameters,
        sinking_tracers,
        open_bottom,
        scalar_type,
        arch,
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
        chlorophyll_ratio=convert(eltype(grid), phytoplankton_chlorophyll_ratio),
    )
end

"""Construct coupled FrankenLOBSTER: Agate living ecology inside OceanBioME LOBSTER."""
function construct(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    parameters::NamedTuple=(;),
    phytoplankton_chlorophyll_ratio=1.31,
    grid=BoxModelGrid(),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    kwargs...,
)
    plankton = _construct_plankton(;
        size_structure, parameters, phytoplankton_chlorophyll_ratio, grid,
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
