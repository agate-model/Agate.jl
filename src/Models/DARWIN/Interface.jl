"""Public constructor for the simplified DARWIN-like ecosystem model.

`DARWIN.construct` builds an Oceananigans/OceanBioME-compatible biogeochemistry
instance for a size-structured plankton model with elemental cycling.

This constructor keeps the surface small:

- structure: choose `n_phyto`, `n_zoo`, and diameter specifications
- parameters: override named parameters via `parameters=(; ...)`
- roles: optionally define consumer/prey membership by group symbols via `roles=(consumers=(...), prey=(...))`
- interaction_overrides: optionally override interaction matrices

For ease of use, interaction overrides are exposed as two separate keywords:

- `palatability_matrix`
- `assimilation_matrix`
"""

using OceanBioME: BoxModelGrid

import ...Utils
import ...Constructor
import ...Interface

export construct

"""
    construct(; kw...) -> bgc

Construct a DARWIN biogeochemistry instance.

Keywords
--------
- `n_phyto=2`, `n_zoo=2`: number of phytoplankton and zooplankton size classes
- `phyto_diameters=(1.5, 20.0, :log_splitting)`: diameter specification for phytoplankton
- `zoo_diameters=(20.0, 100.0, :log_splitting)`: diameter specification for zooplankton
- `parameters=(;)`: parameter overrides (validated against the DARWIN parameter set)
- `roles=nothing`: optional role membership for consumer/prey axes; provide as `(; consumers=(:Z, ...), prey=(:P, ...))`
- `palatability_matrix=nothing`, `assimilation_matrix=nothing`: optional interaction matrices. Each may be:
  - a full `(n_total, n_total)` matrix
  - a rectangular `(n_consumer, n_prey)` matrix (the canonical interaction representation)
  - axis-local group-block matrices sized `(n_consumer_groups, n_prey_groups)`
  - a group-block matrix over *all* groups, wrapped as `GroupBlockMatrix(B)` (expanded during construction)
  - a provider function `(ctx) -> matrix` that returns any of the above
- `grid=BoxModelGrid()`: grid used for precision/architecture inference and sinking velocity fields
- `arch=nothing`: override the architecture (usually inferred from `grid`)
- `sinking_tracers=nothing`: sinking speed overrides, e.g. `(POC = 10/day, ...)`
- `open_bottom=true`: whether sinking tracers leave the domain

Returns
-------
An `Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry` instance.
"""
function construct(;
    n_phyto::Int=2,
    n_zoo::Int=2,
    phyto_diameters=(1.5, 20.0, :log_splitting),
    zoo_diameters=(20.0, 100.0, :log_splitting),
    parameters::NamedTuple=(;),
    roles=nothing,
    palatability_matrix=nothing,
    assimilation_matrix=nothing,
    grid=BoxModelGrid(),
    arch=nothing,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    n_phyto >= 1 || throw(ArgumentError("n_phyto must be >= 1"))
    n_zoo >= 1 || throw(ArgumentError("n_zoo must be >= 1"))

    factory = DarwinFactory()

    base = Interface.default_community(factory)
    community = Constructor.build_plankton_community(
        base;
        n=(Z=n_zoo, P=n_phyto),
        diameters=(
            Z=Utils.diameter_specification(zoo_diameters),
            P=Utils.diameter_specification(phyto_diameters),
        ),
    )

    # Interaction overrides (optional).
    pairs = Pair{Symbol,Any}[]
    palatability_matrix !== nothing &&
        push!(pairs, :palatability_matrix => palatability_matrix)
    assimilation_matrix !== nothing &&
        push!(pairs, :assimilation_matrix => assimilation_matrix)

    interaction_overrides = isempty(pairs) ? nothing : (; pairs...)

    return Constructor.construct_factory(
        factory;
        community=community,
        parameters=parameters,
        roles=roles,
        interaction_overrides=interaction_overrides,
        arch=arch,
        sinking_tracers=sinking_tracers,
        grid=grid,
        open_bottom=open_bottom,
    )
end
