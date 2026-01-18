"""Public constructor for the simplified DARWIN-like ecosystem model.

`DARWIN.construct` builds an Oceananigans/OceanBioME-compatible biogeochemistry
instance for a size-structured plankton model with elemental cycling.

This constructor keeps the surface small:

- structure: choose `n_phyto`, `n_zoo`, and diameter specifications
- dynamics: optionally swap plankton dynamics and (optionally) override selected
  biogeochemical tracer dynamics by key
- parameters: override named parameters via `parameters=(; ...)`
- interactions: optionally override interaction matrices or provide an advanced
  `(ctx) -> NamedTuple` callback
"""

using OceanBioME: BoxModelGrid

using .Tracers:
    phytoplankton_growth_two_nutrients_geider_light,
    zooplankton_default

import ...Utils
import ...Constructor
import ...FactoryInterface

export construct

# Same merge helper as NiPiZD, local to avoid cross-module coupling.
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

Construct a DARWIN biogeochemistry instance.

Keywords
--------
- `n_phyto=2`, `n_zoo=2`: number of phytoplankton and zooplankton size classes
- `phyto_diameters=(1.5, 20.0, :log_splitting)`: diameter specification for phytoplankton
- `zoo_diameters=(20.0, 100.0, :log_splitting)`: diameter specification for zooplankton
- `phyto_dynamics`, `zoo_dynamics`: plankton dynamics builders
- `biogeochem_dynamics=nothing`: optional `NamedTuple` overriding selected tracer dynamics keys
- `parameters=(;)`: parameter overrides (validated against the DARWIN parameter set)
- `palatability_matrix=nothing`, `assimilation_matrix=nothing`: optional interaction matrices
- `interactions=nothing`: optional advanced interaction overrides (NamedTuple)
- `grid=BoxModelGrid()`: grid used for precision/architecture inference and sinking velocity fields
- `arch=nothing`: override the architecture (usually inferred from `grid`)
- `sinking_tracers=nothing`: sinking speed overrides, e.g. `(POC = 10/day, ...)`
- `open_bottom=true`: whether sinking tracers leave the domain

Returns
-------
An `Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry` instance.
"""
function construct(;
    n_phyto::Int = 2,
    n_zoo::Int = 2,
    phyto_diameters = (1.5, 20.0, :log_splitting),
    zoo_diameters = (20.0, 100.0, :log_splitting),
    phyto_dynamics = phytoplankton_growth_two_nutrients_geider_light,
    zoo_dynamics = zooplankton_default,
    biogeochem_dynamics = nothing,
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

    factory = DarwinFactory()
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

    default_bgc = FactoryInterface.default_biogeochem_dynamics(factory)
    merged_bgc = biogeochem_dynamics === nothing ? default_bgc : merge(default_bgc, biogeochem_dynamics)

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
        biogeochem_dynamics = merged_bgc,
        community = community,
        parameters = parameters,
        interactions = merged_interactions,
        arch = arch,
        sinking_tracers = sinking_tracers,
        grid = grid,
        open_bottom = open_bottom,
    )
end
