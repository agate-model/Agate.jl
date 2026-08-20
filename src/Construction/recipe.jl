using ..Factories:
    AbstractBGCFactory,
    DerivedDefault,
    parameter_definitions,
    parameter_directory,
    default_components,
    default_processes
using ..Configuration: normalize_diameters

"""Return the stable recipe-family identifier for a model factory."""
function recipe_family(factory::AbstractBGCFactory)
    throw(ArgumentError("Recipe support is not implemented for $(typeof(factory))."))
end

"""Construct the model factory identified by a recipe family."""
function recipe_factory(::Val{family}) where {family}
    throw(ArgumentError("Unsupported recipe model family $(repr(String(family)))."))
end

"""Authored scientific definition captured before model-state materialization.

`ModelRecipe` stores a stable model-family identifier and the authored scientific
inputs needed to describe a model independently of runtime execution choices such
as grids, architectures, and scalar precision. Model-family code supplies parameter
defaults and tracer equations when the recipe is replayed.
"""
struct ModelRecipe{C,PO,IO,ER,IR,PR,A,S}
    family::Symbol
    community::C
    parameter_overrides::PO
    interaction_overrides::IO
    ecological_roles::ER
    interaction_roles::IR
    parameter_roles::PR
    auxiliary_fields::A
    sinking_tracers::S
    open_bottom::Bool
end

"""Canonical component/process recipe captured before runtime realization.

`ProcessModelRecipe` stores the scientific component/process definition together with
subgroup realization and authored overrides. Runtime precision, architecture, host fields,
parameter materialization, topology maps, and compiled equations are derived on replay.
"""
struct ProcessModelRecipe{C,P,G,CM,PO,IO,S}
    family::Symbol
    components::C
    processes::P
    population_groups::G
    community::CM
    parameter_overrides::PO
    interaction_overrides::IO
    sinking_tracers::S
    open_bottom::Bool
end

"""Resolved deterministic scientific state produced by model construction.

`ModelManifest` records the fully materialized parameters, expanded groups and
tracer ordering, resolved role indices, interaction sources, sinking
configuration, and scalar type. It is an in-memory record of the constructed
model state; durable replay is defined by the corresponding recipe representation.
"""
struct ModelManifest{P,G,TR,A,D,ER,IR,PRI,I,S,T<:Real}
    parameters::P
    group_tracers::G
    tracer_order::TR
    auxiliary_fields::A
    plankton_diameters::D
    ecological_roles::ER
    interaction_role_indices::IR
    parameter_role_indices::PRI
    interaction_matrix_sources::I
    sinking_tracers::S
    open_bottom::Bool
    scalar_type::Type{T}
end

function _structural_isequal(a, b)
    a === b && return true
    typeof(a) === typeof(b) || return false

    if a isa Number || a isa AbstractString || a isa Symbol || a isa Char
        return isequal(a, b)
    elseif a isa AbstractArray
        axes(a) == axes(b) || return false
        return all(_structural_isequal(a[i], b[i]) for i in eachindex(a, b))
    elseif a isa NamedTuple
        keys(a) == keys(b) || return false
        return all(_structural_isequal(getproperty(a, k), getproperty(b, k)) for k in keys(a))
    elseif a isa Tuple
        length(a) == length(b) || return false
        return all(_structural_isequal(a[i], b[i]) for i in eachindex(a))
    elseif isstructtype(typeof(a)) && fieldcount(typeof(a)) > 0
        return all(
            _structural_isequal(getfield(a, i), getfield(b, i))
            for i in 1:fieldcount(typeof(a))
        )
    end

    return isequal(a, b)
end

function _recipe_isequal(a, b)
    a === b && return true

    if a isa Bool || b isa Bool
        return a isa Bool && b isa Bool && isequal(a, b)
    elseif a isa Number || b isa Number
        return a isa Number && b isa Number && isequal(a, b)
    elseif a isa AbstractArray || b isa AbstractArray
        a isa AbstractArray && b isa AbstractArray || return false
        axes(a) == axes(b) || return false
        return all(_recipe_isequal(a[i], b[i]) for i in eachindex(a, b))
    elseif a isa NamedTuple || b isa NamedTuple
        a isa NamedTuple && b isa NamedTuple || return false
        keys(a) == keys(b) || return false
        return all(_recipe_isequal(getproperty(a, k), getproperty(b, k)) for k in keys(a))
    elseif a isa Tuple || b isa Tuple
        a isa Tuple && b isa Tuple || return false
        length(a) == length(b) || return false
        return all(_recipe_isequal(a[i], b[i]) for i in eachindex(a))
    end

    Ta, Tb = typeof(a), typeof(b)
    if isstructtype(Ta) && isstructtype(Tb) && fieldcount(Ta) > 0 &&
       Base.typename(Ta) === Base.typename(Tb) && fieldcount(Ta) == fieldcount(Tb)
        return all(_recipe_isequal(getfield(a, i), getfield(b, i)) for i in 1:fieldcount(Ta))
    end

    return Ta === Tb && isequal(a, b)
end

function Base.:(==)(a::ModelRecipe, b::ModelRecipe)
    return all(
        _recipe_isequal(getfield(a, i), getfield(b, i))
        for i in 1:fieldcount(typeof(a))
    )
end

function Base.:(==)(a::ProcessModelRecipe, b::ProcessModelRecipe)
    return all(
        _recipe_isequal(getfield(a, i), getfield(b, i))
        for i in 1:fieldcount(typeof(a))
    )
end

function Base.:(==)(a::ModelManifest, b::ModelManifest)
    return _structural_isequal(a, b)
end

function _normalize_recipe_community(community::NamedTuple)
    groups = keys(community)
    values = ntuple(length(groups)) do i
        spec = getproperty(community, groups[i])
        diameter_specification = normalize_diameters(spec.diameters).specification
        return (; spec..., diameters=diameter_specification)
    end
    return NamedTuple{groups}(values)
end

"""Capture a construction recipe from semantic inputs before materialization."""
function capture_model_recipe(
    factory::AbstractBGCFactory;
    community,
    parameter_overrides::NamedTuple=(;),
    interaction_overrides::NamedTuple=(;),
    ecological_roles,
    interaction_roles,
    parameter_roles,
    auxiliary_fields::Tuple,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    family = recipe_family(factory)
    family isa Symbol || throw(
        ArgumentError("recipe_family must return a Symbol; got $(typeof(family)).")
    )

    return ModelRecipe(
        family,
        deepcopy(_normalize_recipe_community(community)),
        deepcopy(parameter_overrides),
        deepcopy(interaction_overrides),
        deepcopy(ecological_roles),
        deepcopy(interaction_roles),
        deepcopy(parameter_roles),
        deepcopy(auxiliary_fields),
        deepcopy(sinking_tracers),
        open_bottom,
    )
end

"""Capture a v3 component/process recipe before runtime realization."""
function capture_process_model_recipe(
    factory::AbstractBGCFactory;
    population_groups::NamedTuple,
    community::NamedTuple,
    parameter_overrides::NamedTuple=(;),
    interaction_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    family = recipe_family(factory)
    family isa Symbol || throw(
        ArgumentError("recipe_family must return a Symbol; got $(typeof(family)).")
    )
    return ProcessModelRecipe(
        family,
        deepcopy(default_components(factory)),
        deepcopy(default_processes(factory)),
        deepcopy(population_groups),
        deepcopy(_normalize_recipe_community(community)),
        deepcopy(parameter_overrides),
        deepcopy(interaction_overrides),
        deepcopy(sinking_tracers),
        open_bottom,
    )
end

"""Return the model-family factory used to replay `recipe`."""
replay_factory(recipe::Union{ModelRecipe,ProcessModelRecipe}) = recipe_factory(Val(recipe.family))

"""Capture the resolved deterministic state of a constructed model."""
function capture_model_manifest(
    factory::AbstractBGCFactory,
    parameters,
    community_context;
    tracer_order::Tuple,
    auxiliary_fields::Tuple,
    ecological_roles::NamedTuple=(;),
    explicit_override_keys::Tuple,
    sinking_tracers,
    open_bottom::Bool,
    scalar_type::Type{T},
) where {T<:Real}
    group_order = Tuple(unique(community_context.group_symbols))
    group_values = ntuple(length(group_order)) do i
        group = group_order[i]
        indices = community_context.group_indices[group]
        return Tuple(community_context.plankton_symbols[indices])
    end
    group_tracers = NamedTuple{group_order}(group_values)

    interaction_role_indices = (
        consumers=Tuple(community_context.consumer_indices),
        prey=Tuple(community_context.prey_indices),
    )
    parameter_role_names = keys(community_context.parameter_role_indices)
    parameter_role_values = ntuple(length(parameter_role_names)) do i
        role = parameter_role_names[i]
        Tuple(getproperty(community_context.parameter_role_indices, role))
    end
    parameter_role_indices = NamedTuple{parameter_role_names}(parameter_role_values)

    interaction_names = Tuple(
        spec.name for spec in parameter_directory(factory) if
        spec.shape === :matrix && spec.axes == (:consumer, :prey)
    )
    derived_interaction_names = Tuple(
        definition.spec.name for definition in parameter_definitions(factory) if
        definition.default isa DerivedDefault
    )
    interaction_matrix_sources = NamedTuple{interaction_names}(
        Tuple(
            name in explicit_override_keys ? :explicit :
            name in derived_interaction_names ? :derived : :default
            for name in interaction_names
        ),
    )

    return ModelManifest(
        deepcopy(parameters),
        group_tracers,
        tracer_order,
        auxiliary_fields,
        Tuple(community_context.diameters),
        deepcopy(ecological_roles),
        interaction_role_indices,
        parameter_role_indices,
        deepcopy(interaction_matrix_sources),
        deepcopy(sinking_tracers),
        open_bottom,
        scalar_type,
    )
end
