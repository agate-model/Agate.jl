"""Public constructor for the NiPiZD ecosystem model.

`NiPiZD.construct` builds an Oceananigans/OceanBioME-compatible biogeochemistry
instance for a size-structured plankton model with:

- two plankton groups: phytoplankton (`P`) and zooplankton (`Z`)
- two non-plankton tracers: dissolved inorganic nutrient (`N`) and detritus (`D`)

The public interface keeps the surface small and explicit:

- structure: choose `n_phyto`, `n_zoo`, and diameter specifications
- parameters: override named parameters via `parameters=(; ...)`
- roles: optionally define consumer/prey membership by group symbols via `roles=(consumers=(...), prey=(...))`
- interaction_overrides: optionally override interaction matrices

For ease of use, interaction overrides are exposed as two separate keywords:

- `palatability_matrix`
- `assimilation_matrix`

Each may be either a concrete matrix, or a function that computes a matrix from
the construction context:

- `(community_context) -> matrix`

Matrix overrides may be specified as full `(n_total, n_total)` matrices, or (because
these matrices are role-aware) as rectangular `(n_consumer, n_prey)` matrices.

Internally, role-aware interactions are stored **only** in rectangular form.
No square matrices or square views are created.

To pass a group-block matrix over *all* groups, wrap it as `GroupBlockMatrix(B)` to
force group-block expansion during construction.
Axis-local group-block matrices sized `(n_consumer_groups, n_prey_groups)` are also accepted.
"""

using OceanBioME: BoxModelGrid

import ...Utils
import ...Constructor
import ...Interface

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
- `roles=nothing`: optional role membership for consumer/prey axes; provide as `(; consumers=(:Z, ...), prey=(:P, ...))`
- `parameter_groups=nothing`: optional group membership used only for generating default parameters; provide as `(; producers=(:P, ...), consumers=(:Z, ...))`
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
    roles=nothing,
    parameter_groups=nothing,
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

    roles_resolved = isnothing(roles) ? (consumers=(:Z,), prey=(:P,)) : roles

    parameter_groups_resolved = if isnothing(parameter_groups)
        (
            producers=getproperty(roles_resolved, :prey),
            consumers=getproperty(roles_resolved, :consumers),
        )
    else
        parameter_groups
    end

    return Constructor.construct_factory(
        factory;
        community=community,
        parameters=parameters,
        roles=roles_resolved,
        parameter_groups=parameter_groups_resolved,
        auxiliary_fields=(:PAR,),
        interaction_overrides=interaction_overrides,
        arch=arch,
        sinking_tracers=sinking_tracers,
        grid=grid,
        open_bottom=open_bottom,
    )
end
