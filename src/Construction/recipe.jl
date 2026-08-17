using ..Factories: AbstractBGCFactory
import ..Factories: parameter_definitions
import ..Configuration: matrix_definitions

"""Authored scientific definition captured before model-state materialization.

`ModelRecipe` stores the deterministic construction inputs needed to describe a
model independently of runtime environment objects such as grids and
architectures. Parameter definitions/defaults and user overrides are retained
separately so partial overrides and parameter laws remain semantic inputs rather
than being replaced by resolved vectors.
"""
struct ModelRecipe{
    F,C,PD,PO,ID,IO,GR,IR,DPR,A,S,T<:Real
}
    factory::F
    community::C
    parameter_definitions::PD
    parameter_overrides::PO
    interaction_definitions::ID
    interaction_overrides::IO
    group_roles::GR
    interaction_roles::IR
    default_parameter_roles::DPR
    auxiliary_fields::A
    sinking_tracers::S
    open_bottom::Bool
    scalar_type::Type{T}
end

"""Resolved deterministic model state used for exact construction replay checks."""
struct ModelRealization{P,G,TR,A,D,I,S,T<:Real}
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

function Base.:(==)(a::ModelRecipe, b::ModelRecipe)
    return _structural_isequal(a, b)
end

function Base.:(==)(a::ModelRealization, b::ModelRealization)
    return _structural_isequal(a, b)
end

"""Capture a construction recipe from semantic inputs before materialization."""
function capture_model_recipe(
    factory::AbstractBGCFactory;
    community,
    parameter_overrides::NamedTuple=(;),
    interaction_overrides::Union{Nothing,NamedTuple}=nothing,
    group_roles,
    interaction_roles,
    default_parameter_roles,
    auxiliary_fields::Tuple,
    sinking_tracers=nothing,
    open_bottom::Bool=true,
    scalar_type::Type{T},
) where {T<:Real}
    return ModelRecipe(
        factory,
        deepcopy(community),
        deepcopy(parameter_definitions(factory)),
        deepcopy(parameter_overrides),
        deepcopy(matrix_definitions(factory)),
        deepcopy(interaction_overrides),
        deepcopy(group_roles),
        deepcopy(interaction_roles),
        deepcopy(default_parameter_roles),
        deepcopy(auxiliary_fields),
        deepcopy(sinking_tracers),
        open_bottom,
        scalar_type,
    )
end

"""Internal factory view that reuses the parameter and matrix definitions captured for replay."""
struct ReplayFactory{F,PD,ID} <: AbstractBGCFactory
    factory::F
    parameter_definitions::PD
    interaction_definitions::ID
end

parameter_definitions(factory::ReplayFactory) = factory.parameter_definitions
matrix_definitions(factory::ReplayFactory) = factory.interaction_definitions

replay_factory(recipe::ModelRecipe) = ReplayFactory(
    recipe.factory, recipe.parameter_definitions, recipe.interaction_definitions
)

function capture_model_realization(
    parameters,
    community_context;
    tracer_order::Tuple,
    auxiliary_fields::Tuple,
    interaction_matrix_sources,
    sinking_tracers,
    open_bottom::Bool,
    scalar_type::Type{T},
) where {T<:Real}
    group_order = keys(community_context.plankton_dynamics)
    group_values = ntuple(length(group_order)) do i
        group = group_order[i]
        indices = community_context.group_indices[group]
        return Tuple(community_context.plankton_symbols[indices])
    end
    group_tracers = NamedTuple{group_order}(group_values)

    return ModelRealization(
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
