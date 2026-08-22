using ..ModelFamilies: AbstractModelFamily, default_components, default_processes
using ..Parameters: DerivedDefault, parameter_definitions, parameter_directory
using ..Configuration: normalize_diameters

"""Return the stable recipe-family identifier for a model family."""
function family_id(family::AbstractModelFamily)
    throw(ArgumentError("Recipe support is not implemented for $(typeof(family))."))
end

"""Return the registered model-family token identified by `family`."""
function registered_family(::Val{family}) where {family}
    throw(ArgumentError("Unsupported recipe model family $(repr(String(family)))."))
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
tracer ordering, interaction sources, sinking configuration, and scalar type. It is
an in-memory record of the constructed model state; durable replay is defined by the
corresponding recipe representation.
"""
struct ModelManifest{P,G,TR,A,D,I,S,T<:Real}
    parameters::P
    group_tracers::G
    tracer_order::TR
    auxiliary_fields::A
    plankton_diameters::D
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

"""Capture a component/process recipe before runtime realization."""
function capture_process_model_recipe(
    family::AbstractModelFamily;
    population_groups::NamedTuple,
    community::NamedTuple,
    parameter_overrides::NamedTuple=(;),
    interaction_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    family_id_value = family_id(family)
    family_id_value isa Symbol || throw(
        ArgumentError("family_id must return a Symbol; got $(typeof(family_id_value)).")
    )
    return ProcessModelRecipe(
        family_id_value,
        deepcopy(default_components(family)),
        deepcopy(default_processes(family)),
        deepcopy(population_groups),
        deepcopy(_normalize_recipe_community(community)),
        deepcopy(parameter_overrides),
        deepcopy(interaction_overrides),
        deepcopy(sinking_tracers),
        open_bottom,
    )
end

"""Return the registered model family used to replay `recipe`."""
replay_family(recipe::ProcessModelRecipe) = registered_family(Val(recipe.family))

"""Capture the resolved deterministic state of a constructed model."""
function capture_model_manifest(
    family::AbstractModelFamily,
    parameters,
    community_context;
    tracer_order::Tuple,
    auxiliary_fields::Tuple,
    explicit_override_keys::Tuple,
    sinking_tracers,
    open_bottom::Bool,
    scalar_type::Type{T},
) where {T<:Real}
    group_order = Tuple(unique(community_context.group_symbols))
    group_values = ntuple(length(group_order)) do i
        group = group_order[i]
        indices = community_context.group_indices[group]
        return Tuple(community_context.class_symbols[indices])
    end
    group_tracers = NamedTuple{group_order}(group_values)

    interaction_names = Tuple(
        spec.name for spec in parameter_directory(family) if
        spec.shape === :matrix && spec.axes == (:consumer, :prey)
    )
    derived_interaction_names = Tuple(
        definition.spec.name for definition in parameter_definitions(family) if
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
        deepcopy(interaction_matrix_sources),
        deepcopy(sinking_tracers),
        open_bottom,
        scalar_type,
    )
end
