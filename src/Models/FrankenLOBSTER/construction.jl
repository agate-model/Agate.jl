using OceanBioME: BoxModelGrid
using OceanBioME.Models.NutrientsPlanktonDetritusModels: LOBSTER
using ...Construction

const _SIZE_ROLES = (:phytoplankton, :zooplankton, :bacterioplankton)

function _canonicalize_size_structure(size_structure)
    size_structure isa NamedTuple || throw(
        ArgumentError("size_structure must be a NamedTuple with phytoplankton, zooplankton, and bacterioplankton fields"),
    )
    Set(keys(size_structure)) == Set(_SIZE_ROLES) || throw(
        ArgumentError("size_structure must define exactly phytoplankton, zooplankton, and bacterioplankton"),
    )

    roles = ntuple(length(_SIZE_ROLES)) do i
        role = _SIZE_ROLES[i]
        pfts = getproperty(size_structure, role)
        pfts isa NamedTuple || throw(
            ArgumentError("size_structure.$role must be a NamedTuple of PFT size specifications"),
        )
        isempty(pfts) && throw(
            ArgumentError("size_structure.$role must contain at least one PFT"),
        )
        pfts
    end

    pft_names = Symbol[name for pfts in roles for name in keys(pfts)]
    duplicates = Tuple(name for name in unique(pft_names) if count(==(name), pft_names) > 1)
    isempty(duplicates) || throw(
        ArgumentError("plankton PFT names must be unique across roles; duplicated PFTs: $(collect(duplicates))"),
    )

    return NamedTuple{_SIZE_ROLES}(roles)
end

function _plankton_realization(size_structure)
    structure = _canonicalize_size_structure(size_structure)
    normalize(role) = NamedTuple{keys(getproperty(structure, role))}(Tuple(
        Construction.normalize_pft_size_structure(value)
        for value in values(getproperty(structure, role))
    ))
    return (P=normalize(:phytoplankton), Z=normalize(:zooplankton), H=normalize(:bacterioplankton))
end

"""Construct the internal Agate-generated plankton component consumed by LOBSTER.

The compiled Agate runtime retains OceanBioME-owned resource identities as inputs, while the
returned plankton component registers only realized P/Z/H tracers as living prognostic state.
"""
function _construct_plankton(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    grid=BoxModelGrid(),
    parameters::NamedTuple=(;),
    phytoplankton_chlorophyll_ratio=1.31,
    scalar_type=nothing,
    arch=nothing,
)
    runtime = Construction.construct(
        FrankenLOBSTERFamily();
        plankton_pfts=_plankton_realization(size_structure),
        grid,
        parameter_overrides=parameters,
        scalar_type,
        arch,
    )

    owned = runtime.metadata.plankton_tracers
    phytoplankton = runtime.metadata.pft_entities.P
    chlorophyll_ratio = convert(eltype(grid), phytoplankton_chlorophyll_ratio)
    return FrankenLOBSTERPlankton(
        runtime,
        owned,
        (:solid_waste, :inorganic_waste);
        phytoplankton_tracers=phytoplankton,
        chlorophyll_ratio,
    )
end

"""
    construct(; kw...) -> biogeochemistry

Construct the coupled FrankenLOBSTER model. Agate owns the realized P/Z/H living community
and its biological exchange fluxes; OceanBioME LOBSTER owns nitrate/ammonium, DOM/POM,
remineralization, particle sinking, light, and optional carbon/oxygen components.

`size_structure` and `parameters` configure the Agate living community. Remaining keyword
arguments are forwarded to `OceanBioME.Models.NutrientsPlanktonDetritusModels.LOBSTER`.
"""
function construct(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    parameters::NamedTuple=(;),
    phytoplankton_chlorophyll_ratio=1.31,
    grid=BoxModelGrid(),
    open_bottom::Bool=true,
    kwargs...
)
    plankton = _construct_plankton(;
        size_structure,
        parameters,
        phytoplankton_chlorophyll_ratio,
        grid,
    )

    return LOBSTER(
        grid;
        limiting_nutrients=(:nitrate, :ammonia),
        plankton,
        open_bottom,
        kwargs...
    )
end
