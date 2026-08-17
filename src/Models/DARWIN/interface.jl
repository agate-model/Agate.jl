using OceanBioME: BoxModelGrid

import ...Configuration
import ...Construction
import ...Factories

import ...Manifests
import ...Manifests: construct_from_manifest
using ...Manifests: default_model_manifest

export construct, construct_with_recipe, construct_with_realization, construct_with_manifest

function _community_inputs(phyto_size_structure, zoo_size_structure)
    factory = DarwinFactory()
    base = Factories.default_community(factory)
    community = Configuration.build_plankton_community(
        base; diameters=(Z=zoo_size_structure, P=phyto_size_structure)
    )

    group_order = keys(community)
    recipe_values = ntuple(length(group_order)) do i
        group = group_order[i]
        spec = getproperty(community, group)
        diameter_specification =
            Configuration.normalize_diameters(spec.diameters).specification
        return (; spec..., diameters=diameter_specification)
    end
    recipe_community = Configuration.build_plankton_community(
        NamedTuple{group_order}(recipe_values)
    )

    ecological_roles = (phytoplankton=(:P,), zooplankton=(:Z,))
    interaction_roles = (consumers=(:Z,), prey=(:P,))
    parameter_roles = (producers=(:P,), consumers=(:Z,))

    return (; community, recipe_community, ecological_roles, interaction_roles, parameter_roles)
end

Base.@kwdef struct DARWINConstructionOptions
    phyto_size_structure = (n=2, min_esd=1.5, max_esd=20.0, splitting=:log_splitting)
    zoo_size_structure = (n=2, min_esd=20.0, max_esd=100.0, splitting=:log_splitting)
    parameters::NamedTuple = (;)
    scalar_type = nothing
    palatability_matrix = nothing
    assimilation_matrix = nothing
    grid = BoxModelGrid()
    arch = nothing
    sinking_tracers = nothing
    open_bottom::Bool = true
end

const RECIPE_INPUT_FIELDS = (
    :phyto_size_structure,
    :zoo_size_structure,
    :parameters,
    :scalar_type,
    :palatability_matrix,
    :assimilation_matrix,
    :sinking_tracers,
    :open_bottom,
)
const ENVIRONMENT_INPUT_FIELDS = (:grid, :arch)

_construction_inputs(; kwargs...) = _construction_inputs(DARWINConstructionOptions(; kwargs...))

function _construction_inputs(options::DARWINConstructionOptions)
    (;
        phyto_size_structure,
        zoo_size_structure,
        parameters,
        scalar_type,
        palatability_matrix,
        assimilation_matrix,
        grid,
        arch,
        sinking_tracers,
        open_bottom,
    ) = options

    factory = DarwinFactory()
    community_inputs = _community_inputs(phyto_size_structure, zoo_size_structure)

    pairs = Pair{Symbol,Any}[]
    palatability_matrix !== nothing &&
        push!(pairs, :palatability_matrix => palatability_matrix)
    assimilation_matrix !== nothing &&
        push!(pairs, :assimilation_matrix => assimilation_matrix)

    interaction_overrides = isempty(pairs) ? nothing : (; pairs...)
    resolved_scalar_type = Construction.resolve_construction_scalar_type(grid, scalar_type)
    auxiliary_fields = (:PAR,)

    recipe_kwargs = (;
        community=community_inputs.recipe_community,
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
            plankton_dynamics=default_plankton_dynamics(factory),
            biogeochem_dynamics=default_biogeochem_dynamics(factory),
            community=community_inputs.community,
            parameters=parameters,
            ecological_roles=community_inputs.ecological_roles,
            interaction_roles=community_inputs.interaction_roles,
            parameter_roles=community_inputs.parameter_roles,
            auxiliary_fields,
            interaction_overrides,
            arch,
            sinking_tracers,
            grid,
            scalar_type,
            open_bottom,
        ),
    )
end

"""
    construct(; kw...) -> bgc

Construct a simplified DARWIN-like, size-structured ecosystem model.

!!! info
    This model is in active development and has not been validated against `MITgcm-DARWIN`.

!!! formulation
    TRACERS:

    ∂t cⱼ = ``Uⱼ``DIC - ``Mⱼ`` + ``Gⱼ`` - ``gⱼ``

    ∂t DIC = ∑(``Uⱼ`` DIC) + ``R``DOC + ``R``POC

    ∂t DIN = ∑(``Uⱼ``DIC * ``Qⱼ``N)  + ``R``DON + ``R``PON

    ∂t PO4 = ∑(``Uⱼ``DIC * ``Qⱼ``P)  + ``R``DOP + ``R``POP

    ∂t DOC = ∑(``Mⱼ``DOC) + ``g``DOC - ``R``DOC

    ∂t DON = ∑(``Mⱼ``DOC * ``Qⱼ``N) + ``g``DON - ``R``DON

    ∂t DOP = ∑(``Mⱼ``DOC * ``Qⱼ``P) + ``g``DOP - ``R``DOP

    ∂t POC = ∑(``Mⱼ``POC) + ``g``POC - ``R``POC

    ∂t PON = ∑(``Mⱼ``POC * ``Qⱼ``N) + ``g``PON - ``R``PON

    ∂t POP = ∑(``Mⱼ``POC * ``Qⱼ``P) + ``g``POP - ``R``POP

    where:
    - ``U`` = uptake
    - ``R`` = remineralization
    - ``M`` = mortality
    - ``g, G`` = grazing losses and gains
    - ``Q`` = plankton elemental ratios

    TRAITS:

    μmax, KR, gmax = a*Volume^b

    palat = η/(1+(``ratio``-``opt``)^2)^σ

    where:
    - μmax = maximum photosynthetic growth
    - KR = nutrient half saturation
    - gmax = maximum predation rate
    - palat = palatability
    - ``ratio`` = predator to prey size ratio (diameter)
    - ``opt`` = predator to prey size optimum (diameter)
    - η = prey protection
    - σ = predator specificity

The DARWIN model contains two plankton groups: phytoplankton (`P`) and zooplankton (`Z`),
with size classes defined by their size-structure inputs.

In addition to plankton, the default DARWIN factory includes elemental cycling tracers:
dissolved inorganic carbon (DIC), dissolved inorganic nitrogen (DIN), phosphate (PO4),
dissolved organic matter (DOC, DON, DOP), and particulate organic matter (POC, PON, POP).

During construction, plankton size (diameter) is used to resolve trait-based parameter
vectors and interaction matrices (e.g. palatability and assimilation efficiency). You
may override interaction matrices explicitly with `palatability_matrix` and/or
`assimilation_matrix`.

Size-structure inputs may be a NamedTuple range, for example
`(n=3, min_esd=1, max_esd=10, splitting=:log_splitting)`, or an explicit
diameter vector such as `[1.0, 3.2, 10.0]`.

The returned biogeochemistry instance includes a photosynthetically active radiation (PAR)
auxiliary field.

Keywords
--------
- `phyto_size_structure=(n=2, min_esd=1.5, max_esd=20.0, splitting=:log_splitting)`: phytoplankton size structure
- `zoo_size_structure=(n=2, min_esd=20.0, max_esd=100.0, splitting=:log_splitting)`: zooplankton size structure
- `parameters=(;)`: parameter overrides (validated against the DARWIN parameter set). Vector parameters may be supplied positionally or as partial NamedTuple overrides keyed by plankton class name.
- `scalar_type=nothing`: explicit runtime scalar type. When omitted, construction uses `eltype(grid)` or `Float64` if no grid is supplied
- `palatability_matrix=nothing`: optional palatability matrix override. Must be an explicit rectangular matrix sized to the canonical interaction axes `(n_consumer, n_prey)`.
- `assimilation_matrix=nothing`: optional assimilation matrix override. Must be an explicit rectangular matrix sized to the canonical interaction axes `(n_consumer, n_prey)`

- `grid=BoxModelGrid()`: grid used for precision/architecture inference and sinking velocity fields
- `arch=nothing`: override the architecture (usually inferred from `grid`)
- `sinking_tracers=nothing`: sinking speed overrides, e.g. `(POC = 10/day, ...)`
- `open_bottom=true`: whether sinking tracers leave the domain

Returns
-------
An `Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry` instance.

Example
-------
```julia
using Agate.Models: DARWIN

bgc = DARWIN.construct()
```

Trait-style allometric parameter overrides may be supplied during construction:

```julia
using Oceananigans.Units: day
using Agate.Library.Allometry: AllometricParam, PowerLaw

bgc = DARWIN.construct(;
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

function construct_from_manifest(
    ::Val{:DARWIN}, setup::AbstractDict; grid=nothing, arch=nothing
)
    kwargs = Manifests.manifest_kwargs(
        setup, ("phyto_size_structure", "zoo_size_structure")
    )
    common = Manifests.common_constructor_kwargs(kwargs; grid, arch)
    phyto_size_structure = Manifests.size_structure_vector(kwargs, "phyto_size_structure")
    zoo_size_structure = Manifests.size_structure_vector(kwargs, "zoo_size_structure")
    return construct(; phyto_size_structure, zoo_size_structure, common...)
end

function _recipe_construction_inputs(
    recipe::Construction.ModelRecipe; grid=BoxModelGrid(), arch=nothing
)
    recipe.factory isa DarwinFactory || throw(
        ArgumentError(
            "DARWIN.construct(recipe) requires a DARWIN recipe; got $(typeof(recipe.factory))"
        ),
    )

    factory = Construction.replay_factory(recipe)
    kwargs = (;
        plankton_dynamics=default_plankton_dynamics(recipe.factory),
        biogeochem_dynamics=default_biogeochem_dynamics(recipe.factory),
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

"""Replay a DARWIN recipe in the supplied runtime environment."""
function construct(recipe::Construction.ModelRecipe; grid=BoxModelGrid(), arch=nothing)
    inputs = _recipe_construction_inputs(recipe; grid, arch)
    return Construction.construct_factory(inputs.factory; inputs.kwargs...)
end

"""Construct DARWIN and return the model with its pre-materialization scientific recipe."""
function construct_with_recipe(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    recipe = Construction.capture_model_recipe(inputs.factory; inputs.recipe_kwargs...)
    bgc = Construction.construct_factory(inputs.factory; inputs.kwargs...)
    return bgc, recipe
end

"""
    construct_with_realization(; kw...) -> bgc, recipe, realization
    construct_with_realization(recipe; grid=BoxModelGrid(), arch=nothing) -> bgc, realization

Construct DARWIN together with the authored recipe and deterministic realization used for
exact migration comparison. Passing a recipe replays its captured parameter defaults and
interaction derivations in the supplied runtime environment.
"""
function construct_with_realization(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    recipe = Construction.capture_model_recipe(inputs.factory; inputs.recipe_kwargs...)
    bgc, realization = Construction.construct_factory_with_realization(
        inputs.factory; inputs.kwargs...
    )
    return bgc, recipe, realization
end

function construct_with_realization(
    recipe::Construction.ModelRecipe; grid=BoxModelGrid(), arch=nothing
)
    inputs = _recipe_construction_inputs(recipe; grid, arch)
    bgc, realization = Construction.construct_factory_with_realization(
        inputs.factory; inputs.kwargs...
    )
    return bgc, realization
end

"""
    construct_with_manifest(; kw...) -> bgc, manifest

Construct a model instance and return it with a JSON-compatible manifest of the resolved model setup.
"""
function construct_with_manifest(; kwargs...)
    inputs = _construction_inputs(; kwargs...)
    bgc, manifest_data = Construction.construct_factory_with_manifest_data(
        inputs.factory; inputs.kwargs...
    )
    manifest = default_model_manifest(:DARWIN, manifest_data)
    return bgc, manifest
end
