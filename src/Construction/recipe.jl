using ..ModelFamilies: AbstractModelFamily, definition_version
using ..Parameters: DerivedDefault
using ..Configuration: ModelLayout, canonicalize_diameters

"""Return the stable recipe-family identifier for a model family."""
function family_id(family::AbstractModelFamily)
    throw(ArgumentError("Recipe support is not implemented for $(typeof(family))."))
end

"""Return the registered model-family token identified by `family`."""
function registered_family(::Val{family}) where {family}
    throw(ArgumentError("Unsupported recipe model family $(repr(String(family)))."))
end

"""Versioned registered-family recipe captured before runtime realization.

`ModelRecipe` stores only the registered family identity, its exact scientific
`definition_version`, and canonical construction inputs. Components, processes, parameter
definitions, runtime precision, host fields, and compiled equations are supplied by the
loaded family implementation when the recipe is replayed.
"""
struct ModelRecipe{G,D,P,S}
    family::Symbol
    definition_version::VersionNumber
    plankton_pfts::G
    pft_size_structures::D
    parameter_overrides::P
    sinking_tracers::S
    open_bottom::Bool
end

"""Resolved deterministic scientific state produced by model construction.

`ModelManifest` records the fully materialized parameters, expanded PFT SizeClasses and
tracer ordering, interaction sources, sinking configuration, and scalar type. It is
an in-memory record of the constructed model state; durable replay is defined by the
corresponding recipe representation.
"""
struct ModelManifest{P,G,TR,A,D,I,S,T<:Real}
    parameters::P
    pft_size_classes::G
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

function Base.:(==)(a::ModelRecipe, b::ModelRecipe)
    a.family === b.family || return false
    a.definition_version == b.definition_version || return false
    return _encode_realization(a) == _encode_realization(b)
end

function Base.:(==)(a::ModelManifest, b::ModelManifest)
    return _structural_isequal(a, b)
end

function _canonical_pft_size_structures(pft_size_structures::NamedTuple)
    names = keys(pft_size_structures)
    values = ntuple(length(names)) do i
        pft = names[i]
        canonicalize_diameters(
            getproperty(pft_size_structures, pft); path="plankton PFT :$pft size_structure"
        ).specification
    end
    return NamedTuple{names}(values)
end

"""Capture canonical family construction inputs for durable replay."""
function capture_model_recipe(
    family::AbstractModelFamily;
    plankton_pfts::NamedTuple,
    pft_size_structures::NamedTuple,
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
    return ModelRecipe(
        family_id_value,
        version,
        deepcopy(plankton_pfts),
        deepcopy(_canonical_pft_size_structures(pft_size_structures)),
        deepcopy(parameter_overrides),
        deepcopy(sinking_tracers),
        open_bottom,
    )
end

_family_realization(inputs::NamedTuple) = (;
    plankton_pfts=inputs.plankton_pfts,
    pft_size_structures=inputs.pft_size_structures,
    parameter_overrides=inputs.parameter_overrides,
    sinking_tracers=inputs.sinking_tracers,
    open_bottom=inputs.open_bottom,
)

_family_realization(recipe::ModelRecipe) = (;
    plankton_pfts=recipe.plankton_pfts,
    pft_size_structures=recipe.pft_size_structures,
    parameter_overrides=recipe.parameter_overrides,
    sinking_tracers=recipe.sinking_tracers,
    open_bottom=recipe.open_bottom,
)

_capture_family_recipe(inputs::NamedTuple) =
    capture_model_recipe(inputs.family; _family_realization(inputs)...)

function _resolve_recipe_family(family_id_value::Symbol, version::VersionNumber)
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
replay_family(recipe::ModelRecipe) =
    _resolve_recipe_family(recipe.family, recipe.definition_version)

"""Capture the resolved deterministic state of a constructed model."""
function capture_model_manifest(
    family::AbstractModelFamily,
    parameters,
    layout::ModelLayout,
    parameter_plan;
    interaction_axes,
    tracer_order::Tuple,
    auxiliary_fields::Tuple,
    explicit_override_keys::Tuple,
    sinking_tracers,
    open_bottom::Bool,
    scalar_type::Type{T},
) where {T<:Real}
    pft_order = keys(layout.pft_indices)
    pft_values = ntuple(length(pft_order)) do i
        pft = pft_order[i]
        indices = getproperty(layout.pft_indices, pft)
        return Tuple(layout.size_classes[index] for index in indices)
    end
    pft_size_classes = NamedTuple{pft_order}(pft_values)

    interaction_names = isnothing(interaction_axes) ? () : interaction_axes.parameters
    derived_interaction_names = Tuple(
        name for name in interaction_names if
        getproperty(parameter_plan.parameters, name).definition.default isa DerivedDefault
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
        pft_size_classes,
        tracer_order,
        auxiliary_fields,
        Tuple(layout.size_class_diameters),
        deepcopy(interaction_matrix_sources),
        deepcopy(sinking_tracers),
        open_bottom,
        scalar_type,
    )
end
