export default_components
export default_plankton_dynamics
export default_community
export default_biogeochem_dynamics

"""Canonical logical components for a model factory.

Returns a named collection whose keys are stable model component identities and
whose values describe intrinsic component structure.
"""
function default_components(::AbstractBGCFactory)
    throw(
        ArgumentError(
            "No method `default_components(factory)` is defined for this factory."
        ),
    )
end

"""Default plankton dynamics for a factory.

Returns a `NamedTuple` mapping group symbols (e.g. `:Z`, `:P`) to dynamics builder
functions.
"""
function default_plankton_dynamics(::AbstractBGCFactory)
    throw(
        ArgumentError(
            "No method `default_plankton_dynamics(factory)` is defined for this factory."
        ),
    )
end

"""Resolve plankton dynamics for a realized recipe community.

Factories with fixed group identities can rely on `default_plankton_dynamics(factory)`.
Factories whose authored community uses dynamic group names may specialize this method
using canonical community and ecological-role data.
"""
function default_plankton_dynamics(
    factory::AbstractBGCFactory, community::NamedTuple, ecological_roles::NamedTuple
)
    return default_plankton_dynamics(factory)
end

"""Default plankton community structure for a factory.

Returns a `NamedTuple` mapping group symbols to group specifications.

This is structural information only, such as group symbols, diameter
specifications, and PFT specifications. Numeric parameter defaults are provided
separately through `parameter_definitions(factory)`.
"""
function default_community(::AbstractBGCFactory)
    throw(
        ArgumentError("No method `default_community(factory)` is defined for this factory.")
    )
end

"""Default non-plankton tracer dynamics for a factory.

Returns a `NamedTuple` mapping tracer symbols (e.g. `:N`, `:DIC`) to dynamics
builder functions.
"""
function default_biogeochem_dynamics(::AbstractBGCFactory)
    throw(
        ArgumentError(
            "No method `default_biogeochem_dynamics(factory)` is defined for this factory."
        ),
    )
end
