# """Public constructor for the simplified DARWIN-like ecosystem model.
# 
# `DARWIN.construct` builds an Oceananigans/OceanBioME-compatible biogeochemistry
# instance for a size-structured plankton model with elemental cycling.
# 
# This constructor keeps the surface small:
# 
# - structure: choose `n_phyto`, `n_zoo`, and diameter specifications
# - parameters: override named parameters via `parameters=(; ...)`
# - interaction_overrides: optionally override interaction matrices
# 
# For ease of use, interaction overrides are exposed as two separate keywords:
# 
# - `palatability_matrix`
# - `assimilation_matrix`
# """

using OceanBioME: BoxModelGrid

import ...Configuration
import ...Construction
import ...Factories

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
- `palatability_matrix=nothing`, `assimilation_matrix=nothing`: optional interaction matrix overrides.
  Each must be an explicit rectangular matrix sized to the canonical interaction axes `(n_consumer, n_prey)` (for DARWIN defaults, `(n_zoo, n_phyto)`). Provider functions / callables are not supported; use a Variant/Factory default if you need derived matrices.
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

    base = Factories.default_community(factory)
    community = Configuration.build_plankton_community(
        base;
        n=(Z=n_zoo, P=n_phyto),
        diameters=(
            Z=Configuration.diameter_specification(zoo_diameters),
            P=Configuration.diameter_specification(phyto_diameters),
        ),
    )

    # Interaction overrides (optional).
    # Interaction overrides are data-only: values must be explicit matrices.
    pairs = Pair{Symbol,Any}[]
    palatability_matrix !== nothing &&
        push!(pairs, :palatability_matrix => palatability_matrix)
    assimilation_matrix !== nothing &&
        push!(pairs, :assimilation_matrix => assimilation_matrix)

    interaction_overrides = isempty(pairs) ? nothing : (; pairs...)

    interaction_roles = (consumers=(:Z,), prey=(:P,))

    return Construction.construct_factory(
        factory;
        community=community,
        parameters=parameters,
        interaction_roles=interaction_roles,
        auxiliary_fields=(:PAR,),
        interaction_overrides=interaction_overrides,
        arch=arch,
        sinking_tracers=sinking_tracers,
        grid=grid,
        open_bottom=open_bottom,
    )
end
