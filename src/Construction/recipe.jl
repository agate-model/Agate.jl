using ..ModelFamilies: AbstractModelFamily, default_components
using ..Parameters: DerivedDefault
using ..Components: ModelLayout, canonicalize_plankton_realization, diameter_metadata

"""Return the stable recipe-family identifier for a model family."""
function family_id(family::AbstractModelFamily)
    throw(ArgumentError("Recipe support is not implemented for $(typeof(family))."))
end

"""Canonical construction inputs for an explicitly selected model family.

`ModelRecipe` records the family label and canonical construction configuration.
The selected provider and Julia environment determine the implementation on replay.
Named mapping order is canonicalized for recipe equality and content hashing.
"""
struct ModelRecipe{PlanktonPFTs,ParameterOverrides,SettingOverrides,SinkingTracers}
    family::Symbol
    plankton_pfts::PlanktonPFTs
    parameter_overrides::ParameterOverrides
    setting_overrides::SettingOverrides
    sinking_tracers::SinkingTracers
    open_bottom::Bool
end

"""Resolved deterministic scientific state produced by model construction.

`ModelManifest` records the fully materialized process parameters and model settings, realized
PFT entities and tracer ordering, interaction sources, sinking configuration, and scalar type.
Equality and hashing use this
resolved scientific content; durable replay is defined by the corresponding recipe representation.
"""
struct ModelManifest{
    Parameters,
    Settings,
    PFTEntities,
    TracerOrder,
    AuxiliaryFields,
    PlanktonDiameters,
    InteractionMatrixSources,
    SinkingTracers,
    ScalarType<:Real,
}
    parameters::Parameters
    settings::Settings
    pft_entities::PFTEntities
    tracer_order::TracerOrder
    auxiliary_fields::AuxiliaryFields
    plankton_diameters::PlanktonDiameters
    interaction_matrix_sources::InteractionMatrixSources
    sinking_tracers::SinkingTracers
    open_bottom::Bool
    scalar_type::Type{ScalarType}
end

_recipe_identity(family::Symbol, realization) = (; family, realization)
_recipe_identity(recipe::ModelRecipe) = _recipe_identity(
    recipe.family, _encode_realization(recipe)
)

_manifest_identity(manifest::ModelManifest) = (;
    parameters=manifest.parameters,
    settings=manifest.settings,
    pft_entities=manifest.pft_entities,
    tracer_order=manifest.tracer_order,
    auxiliary_fields=manifest.auxiliary_fields,
    plankton_diameters=manifest.plankton_diameters,
    interaction_matrix_sources=manifest.interaction_matrix_sources,
    sinking_tracers=manifest.sinking_tracers,
    open_bottom=manifest.open_bottom,
    scalar_type=manifest.scalar_type,
)

Base.:(==)(a::ModelRecipe, b::ModelRecipe) = isequal(_recipe_identity(a), _recipe_identity(b))
Base.isequal(a::ModelRecipe, b::ModelRecipe) = isequal(_recipe_identity(a), _recipe_identity(b))
Base.hash(recipe::ModelRecipe, h::UInt) = hash(_recipe_identity(recipe), h)

Base.:(==)(a::ModelManifest, b::ModelManifest) = isequal(_manifest_identity(a), _manifest_identity(b))
Base.isequal(a::ModelManifest, b::ModelManifest) = isequal(_manifest_identity(a), _manifest_identity(b))
Base.hash(manifest::ModelManifest, h::UInt) = hash(_manifest_identity(manifest), h)

function _canonical_recipe_realization(
    family::AbstractModelFamily,
    plankton_pfts::NamedTuple,
)
    return canonicalize_plankton_realization(default_components(family), plankton_pfts)
end


"""Capture canonical family construction inputs for durable replay."""
function capture_model_recipe(
    family::AbstractModelFamily;
    plankton_pfts::NamedTuple,
    parameter_overrides::NamedTuple=(;),
    setting_overrides::NamedTuple=(;),
    sinking_tracers=nothing,
    open_bottom::Bool=true,
)
    family_id_value = family_id(family)
    family_id_value isa Symbol || throw(
        ArgumentError("family_id must return a Symbol; got $(typeof(family_id_value)).")
    )
    plankton_pfts = _canonical_recipe_realization(family, plankton_pfts)
    return ModelRecipe(
        family_id_value,
        deepcopy(plankton_pfts),
        deepcopy(parameter_overrides),
        deepcopy(setting_overrides),
        deepcopy(sinking_tracers),
        open_bottom,
    )
end

_family_realization(recipe::ModelRecipe) = (;
    plankton_pfts=recipe.plankton_pfts,
    parameter_overrides=recipe.parameter_overrides,
    setting_overrides=recipe.setting_overrides,
    sinking_tracers=recipe.sinking_tracers,
    open_bottom=recipe.open_bottom,
)

"""Verify that a recipe belongs to the explicitly selected model family."""
function _check_recipe_family(family::AbstractModelFamily, recipe::ModelRecipe)
    selected = family_id(family)
    selected === recipe.family || throw(ArgumentError(
        "Recipe family $(repr(recipe.family)) does not match selected family $(repr(selected))."
    ))
    return family
end

"""Capture the resolved deterministic state of a constructed model."""
function capture_model_manifest(
    family::AbstractModelFamily,
    parameters,
    settings,
    layout::ModelLayout,
    parameter_plan;
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
    pft_entities = NamedTuple{pft_order}(pft_values)

    interaction_names = Tuple(
        name for (name, parameter) in pairs(parameter_plan.parameters)
        if parameter.runtime_bound && parameter.axes == (:consumer, :resource)
    )
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
        deepcopy(settings),
        pft_entities,
        tracer_order,
        auxiliary_fields,
        Tuple(diameter_metadata(diameter) for diameter in layout.size_class_diameters),
        deepcopy(interaction_matrix_sources),
        deepcopy(sinking_tracers),
        open_bottom,
        scalar_type,
    )
end
