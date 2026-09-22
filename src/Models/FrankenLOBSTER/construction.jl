using OceanBioME: BoxModelGrid
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers

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
    return (P=normalize(:phytoplankton), Z=normalize(:zooplankton), B=normalize(:bacterioplankton))
end

"""Construct the internal Agate-generated plankton component consumed by LOBSTER.

The compiled Agate runtime retains OceanBioME-owned resource identities as inputs, while the
returned plankton component registers only realized P/Z/B tracers as living prognostic state.
"""
function _construct_plankton(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    grid=BoxModelGrid(),
    scalar_type=nothing,
    arch=nothing,
)
    runtime = Construction.construct(
        FrankenLOBSTERFamily();
        plankton_pfts=_plankton_realization(size_structure),
        grid,
        scalar_type,
        arch,
    )

    all_tracers = required_biogeochemical_tracers(runtime)
    owned = runtime.metadata.plankton_tracers
    external = Tuple(tracer for tracer in all_tracers if !(tracer in owned))

    return FrankenLOBSTERPlankton(runtime, owned, external)
end
