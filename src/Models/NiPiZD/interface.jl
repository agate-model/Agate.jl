using OceanBioME: BoxModelGrid

import ...Configuration
import ...Construction

export construct, construct_plus_recipe, construct_from_recipe

function _validated_size_structure(size_structure)
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
        throw(ArgumentError("size_structure.phytoplankton must define at least one group"))
    isempty(zooplankton) &&
        throw(ArgumentError("size_structure.zooplankton must define at least one group"))

    producer_groups = keys(phytoplankton)
    consumer_groups = keys(zooplankton)
    duplicate_groups = [group for group in producer_groups if group in consumer_groups]
    isempty(duplicate_groups) || throw(
        ArgumentError(
            "plankton group names must be unique across roles; " *
            "duplicated groups: $(collect(duplicate_groups))",
        ),
    )

    return (; phytoplankton, zooplankton, producer_groups, consumer_groups)
end

function _community_inputs(size_structure)
    structure = _validated_size_structure(size_structure)
    group_order = (structure.consumer_groups..., structure.producer_groups...)
    empty_pft = Configuration.PFTSpecification()

    community_base_values = ntuple(length(group_order)) do i
        group = group_order[i]
        diameters = if group in structure.consumer_groups
            getproperty(structure.zooplankton, group)
        else
            getproperty(structure.phytoplankton, group)
        end
        return (; diameters, pft=empty_pft)
    end
    community_base = NamedTuple{group_order}(community_base_values)
    community = Configuration.build_plankton_community(community_base)

    dynamics_values = ntuple(length(group_order)) do i
        group = group_order[i]
        return group in structure.consumer_groups ? zooplankton_nipizd : phytoplankton_nipizd
    end
    plankton_dynamics = NamedTuple{group_order}(dynamics_values)
    interaction_roles = (consumers=structure.consumer_groups, prey=structure.producer_groups)
    parameter_roles = (;
        producers=structure.producer_groups, consumers=structure.consumer_groups
    )

    ecological_roles = (;
        phytoplankton=structure.producer_groups, zooplankton=structure.consumer_groups
    )
    return (;
        community,
        plankton_dynamics,
        interaction_roles,
        parameter_roles,
        ecological_roles,
    )
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

    factory = NiPiZDFactory()
    community_inputs = _community_inputs(size_structure)

    pairs = Pair{Symbol,Any}[]
    palatability_matrix !== nothing &&
        push!(pairs, :palatability_matrix => palatability_matrix)
    assimilation_matrix !== nothing &&
        push!(pairs, :assimilation_matrix => assimilation_matrix)

    interaction_overrides = (; pairs...)
    resolved_scalar_type = Construction.resolve_construction_scalar_type(grid, scalar_type)
    auxiliary_fields = (:PAR,)

    recipe_kwargs = (;
        community=community_inputs.community,
        parameter_overrides=parameters,
        interaction_overrides,
        ecological_roles=community_inputs.ecological_roles,
        interaction_roles=community_inputs.interaction_roles,
        parameter_roles=community_inputs.parameter_roles,
        auxiliary_fields,
        sinking_tracers,
        open_bottom,
        scalar_type=resolved_scalar_type,
    )

    return (
        factory=factory,
        recipe_kwargs,
        kwargs=(;
            plankton_dynamics=community_inputs.plankton_dynamics,
            community=community_inputs.community,
            parameters=parameters,
            ecological_roles=community_inputs.ecological_roles,
            interaction_roles=community_inputs.interaction_roles,
            parameter_roles=community_inputs.parameter_roles,
            auxiliary_fields,
            interaction_overrides=interaction_overrides,
            arch=arch,
            sinking_tracers=sinking_tracers,
            grid=grid,
            scalar_type=scalar_type,
            open_bottom=open_bottom,
        ),
    )
end

"""
    construct(; kw...) -> bgc

Construct a size-structured NiPiZD ecosystem model.

The NiPiZD model contains phytoplankton and zooplankton roles. `size_structure` defines
the groups within each role; the defaults are `P` and `Z`.

In addition to plankton, the default NiPiZD factory includes idealized nutrient (`N`) and
detritus (`D`) cycling. The returned biogeochemistry instance includes a photosynthetically
active radiation (PAR) auxiliary field.

During construction, plankton size (diameter) is used to resolve trait-based parameter
vectors and interaction matrices (e.g. palatability and assimilation efficiency). You
may override interaction matrices explicitly with `palatability_matrix` and/or
`assimilation_matrix`.

Each group size structure may be a NamedTuple range, for example
`(n=3, min_esd=1, max_esd=10, splitting=:log_splitting)`, or an explicit
diameter vector such as `[1.0, 3.2, 10.0]`. Groups are supplied as
`size_structure=(phytoplankton=(...), zooplankton=(...))` and classes use
`<group>_<index>` tracer names, such as `P_1` or `diat_1`.

Keywords
--------
- `size_structure`: phytoplankton and zooplankton groups, supplied as a NamedTuple with
  `phytoplankton` and `zooplankton` fields. Defaults to `P` and `Z` groups with two size classes
  each.
- `parameters=(;)`: parameter overrides (validated against the NiPiZD parameter set). Vector parameters may be supplied positionally, as partial NamedTuple overrides keyed by realized plankton tracer name (for example `P_1`, `diat_1`, or `microzoo_1`), or as allometric definitions for diameter-indexed plankton vectors.
- `palatability_matrix=nothing`: optional palatability matrix override. Must be an explicit rectangular matrix with rows ordered by realized zooplankton classes and columns ordered by realized phytoplankton classes.
- `assimilation_matrix=nothing`: optional assimilation matrix override with the same consumer-by-prey class ordering as `palatability_matrix`.
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
    inputs = _construction_inputs(; kwargs...)
    return Construction.construct_factory(inputs.factory; inputs.kwargs...)
end

function _recipe_plankton_dynamics(recipe::Construction.ModelRecipe)
    groups = keys(recipe.community)
    phytoplankton = recipe.ecological_roles.phytoplankton
    zooplankton = recipe.ecological_roles.zooplankton

    values = ntuple(length(groups)) do i
        group = groups[i]
        if group in phytoplankton
            phytoplankton_nipizd
        elseif group in zooplankton
            zooplankton_nipizd
        else
            throw(ArgumentError("recipe group :$group has no NiPiZD ecological role"))
        end
    end
    return NamedTuple{groups}(values)
end

function _recipe_construction_inputs(
    recipe::Construction.ModelRecipe; grid=BoxModelGrid(), arch=nothing
)
    family = recipe.family
    family == :NiPiZD || throw(
        ArgumentError(
            "NiPiZD.construct_from_recipe requires a NiPiZD recipe; got family $family"
        ),
    )

    factory = Construction.replay_factory(recipe)
    kwargs = (;
        plankton_dynamics=_recipe_plankton_dynamics(recipe),
        biogeochem_dynamics=default_biogeochem_dynamics(factory),
        community=recipe.community,
        parameters=recipe.parameter_overrides,
        interaction_overrides=recipe.interaction_overrides,
        ecological_roles=recipe.ecological_roles,
        interaction_roles=recipe.interaction_roles,
        parameter_roles=recipe.parameter_roles,
        auxiliary_fields=recipe.auxiliary_fields,
        arch,
        sinking_tracers=recipe.sinking_tracers,
        grid,
        scalar_type=recipe.scalar_type,
        open_bottom=recipe.open_bottom,
    )
    return (; factory, kwargs)
end

"""
    construct_from_recipe(recipe; grid=BoxModelGrid(), arch=nothing) -> bgc

Replay a NiPiZD recipe in the supplied runtime environment.
"""
function construct_from_recipe(
    recipe::Construction.ModelRecipe; grid=BoxModelGrid(), arch=nothing
)
    inputs = _recipe_construction_inputs(recipe; grid, arch)
    return Construction.construct_factory(inputs.factory; inputs.kwargs...)
end


"""
    construct_plus_recipe(; kw...) -> bgc, recipe

Construct NiPiZD and return the model together with its authored scientific recipe.
The recipe records semantic size specifications, authored parameter and interaction
overrides, role selections, sinking configuration, open-bottom state, and resolved
scalar type. Model-family code supplies defaults, derivations, and equations when the
recipe is replayed. Runtime grid and architecture objects remain construction
environment inputs and are not stored in the recipe.
"""
function construct_plus_recipe(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    recipe = Construction.capture_model_recipe(inputs.factory; inputs.recipe_kwargs...)
    bgc = Construction.construct_factory(inputs.factory; inputs.kwargs...)
    return bgc, recipe
end
