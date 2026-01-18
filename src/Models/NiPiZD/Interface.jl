"""Public constructor for the NiPiZD ecosystem model.

`NiPiZD.construct` builds an Oceananigans/OceanBioME-compatible biogeochemistry
instance for a size-structured plankton model with:

- two plankton groups: phytoplankton (`P`) and zooplankton (`Z`)
- two non-plankton tracers: dissolved inorganic nutrient (`N`) and detritus (`D`)

The public interface keeps the surface small and explicit:

- structure: choose `n_phyto`, `n_zoo`, and diameter specifications
- dynamics: optionally swap any of the four dynamics builders
- parameters: override named parameters via `parameters=(; ...)`
- interactions: optionally override interaction matrices

Advanced users can still pass an `interactions` NamedTuple or a callback
`(ctx) -> NamedTuple`, but the default workflow is to use explicit keyword overrides.
"""

using OceanBioME: BoxModelGrid

using .Tracers:
    nutrient_default,
    detritus_default,
    phytoplankton_default,
    zooplankton_default

import ...Utils
import ...Constructor
import ...FactoryInterface

export construct

# Merge interaction overrides while giving explicit matrices the highest priority.
function _merge_interactions(interactions, matrix_overrides::NamedTuple)
    isempty(keys(matrix_overrides)) && return interactions

    if interactions === nothing
        return matrix_overrides
    elseif interactions isa NamedTuple
        return merge(interactions, matrix_overrides)
    elseif interactions isa Function
        return (ctx) -> merge(interactions(ctx), matrix_overrides)
    else
        throw(ArgumentError("interactions must be nothing, a NamedTuple, or a function (ctx)->NamedTuple"))
    end
end

"""
    construct(; kw...) -> bgc

Construct a NiPiZD biogeochemistry instance.

Keywords
--------
- `n_phyto=2`, `n_zoo=2`: number of phytoplankton and zooplankton size classes
- `phyto_diameters=(2, 10, :log_splitting)`: diameter specification for phytoplankton
- `zoo_diameters=(20, 100, :linear_splitting)`: diameter specification for zooplankton
- `nutrient_dynamics`, `detritus_dynamics`, `phyto_dynamics`, `zoo_dynamics`: dynamics builders
- `parameters=(;)`: parameter overrides (validated against the NiPiZD parameter set)
- `palatability_matrix=nothing`, `assimilation_matrix=nothing`: optional interaction matrices
- `interactions=nothing`: optional advanced interaction overrides (NamedTuple)
- `grid=BoxModelGrid()`: grid used for precision/architecture inference and sinking velocity fields
- `arch=nothing`: override the architecture (usually inferred from `grid`)
- `sinking_tracers=nothing`: sinking speed overrides, e.g. `(D = 2/ day, P1 = 0.1/day, ...)`
- `open_bottom=true`: whether sinking tracers leave the domain

Returns
-------
An `Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry` instance.
"""
function construct(;
    n_phyto::Int = 2,
    n_zoo::Int = 2,
    phyto_diameters = (2, 10, :log_splitting),
    zoo_diameters = (20, 100, :linear_splitting),
    nutrient_dynamics = nutrient_default,
    detritus_dynamics = detritus_default,
    phyto_dynamics = phytoplankton_default,
    zoo_dynamics = zooplankton_default,
    parameters::NamedTuple = (;),
    palatability_matrix = nothing,
    assimilation_matrix = nothing,
    interactions = nothing,
    grid = BoxModelGrid(),
    arch = nothing,
    sinking_tracers = nothing,
    open_bottom::Bool = true,
)
    n_phyto >= 1 || throw(ArgumentError("n_phyto must be >= 1"))
    n_zoo >= 1 || throw(ArgumentError("n_zoo must be >= 1"))

    factory = NiPiZDFactory()
    spec = Constructor.ModelSpec(factory)

    base = FactoryInterface.default_community(factory)
    community = Constructor.build_ZP_community(
        base;
        n_zoo = n_zoo,
        n_phyto = n_phyto,
        zoo_diameters = Utils.diameter_specification(zoo_diameters),
        phyto_diameters = Utils.diameter_specification(phyto_diameters),
    )

    plankton_dynamics = (Z = zoo_dynamics, P = phyto_dynamics)
    biogeochem_dynamics = (N = nutrient_dynamics, D = detritus_dynamics)

    matrix_overrides = NamedTuple()
    if palatability_matrix !== nothing
        matrix_overrides = merge(matrix_overrides, (; palatability_matrix = palatability_matrix))
    end
    if assimilation_matrix !== nothing
        matrix_overrides = merge(matrix_overrides, (; assimilation_matrix = assimilation_matrix))
    end

    merged_interactions = _merge_interactions(interactions, matrix_overrides)

    return Constructor.construct(
        spec;
        plankton_dynamics = plankton_dynamics,
        biogeochem_dynamics = biogeochem_dynamics,
        community = community,
        parameters = parameters,
        interactions = merged_interactions,
        arch = arch,
        sinking_tracers = sinking_tracers,
        grid = grid,
        open_bottom = open_bottom,
    )
end
