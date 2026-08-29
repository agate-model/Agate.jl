using OceanBioME: BoxModelGrid

import ...Configuration
import ...Construction


function _canonicalize_size_structure(size_structure)
    size_structure isa NamedTuple ||
        throw(ArgumentError("size_structure must be a NamedTuple"))

    required_roles = (:phytoplankton, :zooplankton)
    missing_roles = [role for role in required_roles if !hasproperty(size_structure, role)]
    extra_roles = [role for role in keys(size_structure) if !(role in required_roles)]
    isempty(missing_roles) || throw(
        ArgumentError("size_structure is missing roles: $(collect(missing_roles))")
    )
    isempty(extra_roles) ||
        throw(ArgumentError("size_structure has unknown roles: $(collect(extra_roles))"))

    phytoplankton = size_structure.phytoplankton
    zooplankton = size_structure.zooplankton
    phytoplankton isa NamedTuple ||
        throw(ArgumentError("size_structure.phytoplankton must be a NamedTuple"))
    zooplankton isa NamedTuple ||
        throw(ArgumentError("size_structure.zooplankton must be a NamedTuple"))
    isempty(phytoplankton) &&
        throw(ArgumentError("size_structure.phytoplankton must define at least one PFT"))
    isempty(zooplankton) &&
        throw(ArgumentError("size_structure.zooplankton must define at least one PFT"))

    producer_pfts = keys(phytoplankton)
    consumer_pfts = keys(zooplankton)
    duplicate_pfts = [pft for pft in producer_pfts if pft in consumer_pfts]
    isempty(duplicate_pfts) || throw(
        ArgumentError(
            "plankton PFT names must be unique across roles; " *
            "duplicated PFTs: $(collect(duplicate_pfts))",
        ),
    )

    return (; phytoplankton, zooplankton, producer_pfts, consumer_pfts)
end

function _plankton_realization(size_structure)
    structure = _canonicalize_size_structure(size_structure)
    pft_order = (structure.producer_pfts..., structure.consumer_pfts...)
    pft_size_structures = NamedTuple{pft_order}(ntuple(length(pft_order)) do i
        pft = pft_order[i]
        diameters = if pft in structure.consumer_pfts
            getproperty(structure.zooplankton, pft)
        else
            getproperty(structure.phytoplankton, pft)
        end
        diameters
    end)
    plankton_pfts = (P=structure.producer_pfts, Z=structure.consumer_pfts)
    return (; plankton_pfts, pft_size_structures)
end

function _construction_inputs(;
    size_structure=DEFAULT_SIZE_STRUCTURE,
    parameters::NamedTuple=(;),
    palatability_matrix=nothing,
    assimilation_matrix=nothing,
    grid=BoxModelGrid(),
    scalar_type=nothing,
    arch=nothing,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    family = NiPiZDFamily()
    realization = _plankton_realization(size_structure)

    parameter_overrides = parameters
    palatability_matrix === nothing ||
        (parameter_overrides = merge(parameter_overrides, (; palatability_matrix)))
    assimilation_matrix === nothing ||
        (parameter_overrides = merge(parameter_overrides, (; assimilation_matrix)))
    return (;
        family,
        plankton_pfts=realization.plankton_pfts,
        pft_size_structures=realization.pft_size_structures,
        parameter_overrides,
        sinking_tracers,
        open_bottom,
        execution=(; grid, scalar_type, arch),
    )
end

"""
    construct(; kw...) -> bgc

Construct a size-structured NiPiZD ecosystem model.

The NiPiZD model contains phytoplankton and zooplankton roles. `size_structure` defines
the PFTs within each role; the defaults are `P` and `Z`.

In addition to plankton, the NiPiZD definition includes idealized nutrient (`N`) and
detritus (`D`) cycling. The returned biogeochemistry instance includes a photosynthetically
active radiation (PAR) auxiliary field.

During construction, plankton size (diameter) is used to resolve trait-based parameter
vectors and interaction matrices (e.g. palatability and assimilation efficiency). You
may override interaction matrices explicitly with `palatability_matrix` and/or
`assimilation_matrix`.

Each PFT size structure may be a NamedTuple range, for example
`(n=3, min_esd=1, max_esd=10, splitting=:log_splitting)`, or an explicit
diameter vector such as `[1.0, 3.2, 10.0]`. PFTs are supplied as
`size_structure=(phytoplankton=(...), zooplankton=(...))` and SizeClasses use
`<pft>_<index>` identities, such as `P_1` or `diat_1`.

Keywords
--------
- `size_structure`: phytoplankton and zooplankton PFTs, supplied as a NamedTuple with
  `phytoplankton` and `zooplankton` fields. Defaults to `P` and `Z` PFTs with two SizeClasses
  each.
- `parameters=(;)`: parameter overrides (validated against the NiPiZD parameter set). Vector parameters may be supplied positionally, as partial NamedTuple overrides keyed by realized plankton SizeClass identity (for example `P_1`, `diat_1`, or `microzoo_1`), or as allometric definitions for diameter-indexed plankton vectors.
- `palatability_matrix=nothing`: optional palatability matrix override. Must be an explicit rectangular matrix with rows ordered by realized zooplankton SizeClasses and columns ordered by realized phytoplankton SizeClasses.
- `assimilation_matrix=nothing`: optional assimilation matrix override with the same consumer-by-prey SizeClass ordering as `palatability_matrix`.
- `grid=BoxModelGrid()`: grid used for architecture inference and default scalar-type selection
- `scalar_type=nothing`: explicit runtime scalar type. When omitted, construction uses `eltype(grid)` or `Float64` if no grid is supplied
- `arch=nothing`: override the architecture (usually inferred from `grid`)
- `sinking_tracers=nothing`: sinking speed overrides, e.g. `(D = 2/day, P_1 = 0.1/day, ...)`
- `open_bottom=true`: whether sinking tracers leave the domain

Returns
-------
An `Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry` instance.

Example
-------
```julia
using Agate.Models: NiPiZD

bgc = NiPiZD.construct()
```

Trait-style allometric parameter overrides may be supplied during construction:

```julia
using Oceananigans.Units: day
using Agate.Library.Allometry: AllometricParam, PowerLaw

bgc = NiPiZD.construct(;
    parameters=(;
        maximum_growth_rate=AllometricParam(
            PowerLaw(); prefactor=2 / day, exponent=-0.15
        ),
    ),
)
```
"""
function construct(; kwargs...)
    return Construction._construct_family(_construction_inputs(; kwargs...))
end

"""
    construct_from_recipe(recipe; grid=BoxModelGrid(), arch=nothing, scalar_type=nothing) -> bgc

Replay a NiPiZD recipe in the supplied runtime environment.
"""
function construct_from_recipe(
    recipe::Construction.ModelRecipe; grid=BoxModelGrid(), arch=nothing, scalar_type=nothing
)
    recipe.family == :NiPiZD || throw(
        ArgumentError(
            "NiPiZD.construct_from_recipe requires a NiPiZD recipe; got family $(recipe.family)"
        ),
    )
    return Construction.construct(recipe; grid, arch, scalar_type)
end


"""
    construct_plus_recipe(; kw...) -> bgc, recipe

Construct NiPiZD and return the model together with its versioned family recipe.
The recipe records the registered family version, PFT size realization, authored
parameter overrides, sinking configuration, and open-bottom state. Runtime grid,
architecture, scalar precision, components, processes, and compiled equations are supplied
when the recipe is realized.
"""
function construct_plus_recipe(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    recipe = Construction._capture_family_recipe(inputs)
    return Construction._construct_family(inputs), recipe
end
