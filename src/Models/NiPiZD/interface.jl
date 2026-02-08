# """Public constructor for the NiPiZD ecosystem model.
# 
# `NiPiZD.construct` builds an Oceananigans/OceanBioME-compatible biogeochemistry
# instance for a size-structured plankton model with:
# 
# - two plankton groups: phytoplankton (`P`) and zooplankton (`Z`)
# - two non-plankton tracers: dissolved inorganic nutrient (`N`) and detritus (`D`)
# 
# The public interface keeps the surface small and explicit:
# 
# - structure: choose `n_phyto`, `n_zoo`, and diameter specifications
# - parameters: override named parameters via `parameters=(; ...)`
# - interaction_overrides: optionally override interaction matrices
# 
# For ease of use, interaction overrides are exposed as two separate keywords:
# 
# - `palatability_matrix`
# - `assimilation_matrix`
# 
# Each may be either a concrete matrix, or a function that computes a matrix from
# the construction context:
# 
# - `(community_context) -> matrix`
# 
# Matrix overrides may be specified as full `(n_total, n_total)` matrices, or (because
# these matrices are role-aware) as rectangular `(n_consumer, n_prey)` matrices.
# 
# Internally, role-aware interactions are stored **only** in rectangular form.
# No square matrices or square views are created.
# 
# To pass a group-block matrix over *all* groups, wrap it as `GroupBlockMatrix(B)` to
# force group-block expansion during construction.
# Axis-local group-block matrices sized `(n_consumer_groups, n_prey_groups)` are also accepted.
# """

using OceanBioME: BoxModelGrid

import ...Configuration
import ...Construction
import ...Factories

export construct

"""
    construct(; kw...) -> bgc

Construct a NiPiZD biogeochemistry instance.

Keywords
--------
- `n_phyto=2`, `n_zoo=2`: number of phytoplankton and zooplankton size classes
- `phyto_diameters=(2, 10, :log_splitting)`: diameter specification for phytoplankton
- `zoo_diameters=(20, 100, :linear_splitting)`: diameter specification for zooplankton
- `parameters=(;)`: parameter overrides (validated against the NiPiZD parameter set)
- `palatability_matrix=nothing`, `assimilation_matrix=nothing`: optional interaction matrix overrides
- `grid=BoxModelGrid()`: grid used for precision/architecture inference and sinking velocity fields
- `arch=nothing`: override the architecture (usually inferred from `grid`)
- `sinking_tracers=nothing`: sinking speed overrides, e.g. `(D = 2/ day, P1 = 0.1/day, ...)`
- `open_bottom=true`: whether sinking tracers leave the domain

Returns
-------
An `Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry` instance.
"""
function construct(;
    n_phyto::Int=2,
    n_zoo::Int=2,
    phyto_diameters=(2, 10, :log_splitting),
    zoo_diameters=(20, 100, :linear_splitting),
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

    factory = NiPiZDFactory()

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
    #
    # We forward overrides through the model-agnostic constructor as a `NamedTuple`.
    # If a value is a function, it will be evaluated once against the CommunityContext
    # during construction.
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
