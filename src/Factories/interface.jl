export default_plankton_dynamics
export default_community
export default_biogeochem_dynamics

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
