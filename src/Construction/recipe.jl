using ..ModelFamilies: AbstractModelFamily, definition_version
using ..Parameters: DerivedDefault, parameter_definitions, parameter_directory
using ..Configuration:
    PFTSpecification, build_plankton_community, normalize_diameters

"""Return the stable recipe-family identifier for a model family."""
function family_id(family::AbstractModelFamily)
    throw(ArgumentError("Recipe support is not implemented for $(typeof(family))."))
end

"""Return the registered model-family token identified by `family`."""
function registered_family(::Val{family}) where {family}
    throw(ArgumentError("Unsupported recipe model family $(repr(String(family)))."))
end

"""Versioned registered-family recipe captured before runtime realization.

`ProcessModelRecipe` stores only the registered family identity, its exact scientific
`definition_version`, and canonical construction inputs. Components, processes, parameter
definitions, runtime precision, host fields, and compiled equations are supplied by the
loaded family implementation when the recipe is replayed.
"""
struct ProcessModelRecipe{G,D,P,S}
    family::Symbol
    definition_version::VersionNumber
    population_groups::G
    group_diameters::D
    parameter_overrides::P
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

function Base.:(==)(a::ProcessModelRecipe, b::ProcessModelRecipe)
    a.family === b.family || return false
    a.definition_version == b.definition_version || return false
    return _encode_realization(a) == _encode_realization(b)
end

function Base.:(==)(a::ModelManifest, b::ModelManifest)
    return _structural_isequal(a, b)
end

function _canonical_group_diameters(community::NamedTuple)
    names = keys(community)
    values = ntuple(length(names)) do i
        spec = getproperty(community, names[i])
        normalize_diameters(spec.diameters).specification
    end
    return NamedTuple{names}(values)
end

function _recipe_community(recipe::ProcessModelRecipe)
    names = keys(recipe.group_diameters)
    pft = PFTSpecification()
    base = NamedTuple{names}(ntuple(length(names)) do i
        (; diameters=getproperty(recipe.group_diameters, names[i]), pft)
    end)
    return build_plankton_community(base)
end

"""Capture canonical family construction inputs for durable replay."""
function capture_process_model_recipe(
    family::AbstractModelFamily;
    population_groups::NamedTuple,
    community::NamedTuple,
    parameter_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    family_id_value = family_id(family)
    family_id_value isa Symbol || throw(
        ArgumentError("family_id must return a Symbol; got $(typeof(family_id_value)).")
    )
    version = definition_version(family)
    version isa VersionNumber || throw(
        ArgumentError("definition_version must return a VersionNumber; got $(typeof(version)).")
    )
    return ProcessModelRecipe(
        family_id_value,
        version,
        deepcopy(population_groups),
        deepcopy(_canonical_group_diameters(community)),
        deepcopy(parameter_overrides),
        deepcopy(sinking_tracers),
        open_bottom,
    )
end

_family_realization(inputs::NamedTuple) = (;
    population_groups=inputs.population_groups,
    community=inputs.community,
    parameter_overrides=inputs.parameter_overrides,
    sinking_tracers=inputs.sinking_tracers,
    open_bottom=inputs.open_bottom,
)

_family_realization(recipe::ProcessModelRecipe) = (;
    population_groups=recipe.population_groups,
    community=_recipe_community(recipe),
    parameter_overrides=recipe.parameter_overrides,
    sinking_tracers=recipe.sinking_tracers,
    open_bottom=recipe.open_bottom,
)

_capture_family_recipe(inputs::NamedTuple) =
    capture_process_model_recipe(inputs.family; _family_realization(inputs)...)

function _validated_recipe_family(family_id_value::Symbol, version::VersionNumber)
    family = registered_family(Val(family_id_value))
    loaded_family_id = family_id(family)
    loaded_family_id === family_id_value || throw(
        ArgumentError(
            "Registered recipe family $(repr(family_id_value)) resolves to family id " *
            "$(repr(loaded_family_id))."
        )
    )
    loaded_version = definition_version(family)
    version == loaded_version || throw(
        ArgumentError(
            "Recipe definition_version $version does not match loaded " *
            "$family_id_value definition_version $loaded_version."
        )
    )
    return family
end

"""Return the registered model family used to replay `recipe` after identity/version checks."""
replay_family(recipe::ProcessModelRecipe) =
    _validated_recipe_family(recipe.family, recipe.definition_version)

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
        name for (name, spec) in pairs(parameter_directory(family)) if
        spec.axes == (:consumer, :prey)
    )
    derived_interaction_names = Tuple(
        name for (name, parameter) in pairs(parameter_definitions(family)) if
        parameter.default isa DerivedDefault
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
