using ..Factories: AbstractBGCFactory, parameter_definitions
using ..Configuration: matrix_definitions

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
