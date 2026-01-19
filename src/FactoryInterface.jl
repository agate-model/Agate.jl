"""Factory-facing interface for model defaults.

Agate's public model constructors (`NiPiZD.construct`, `DARWIN.construct`, …) are
implemented as thin wrappers around factories. Factories provide *explicit*
defaults for community structure, tracer dynamics, and group ordering.

These functions intentionally live outside `Agate.Models` to avoid circular load
dependencies between the generic constructor pipeline and model modules.
"""

module FactoryInterface

using ..Utils:
    AbstractBGCFactory, ParameterSpec, parameter_directory, consumer_groups, prey_groups

export ParameterSpec
export parameter_directory

export factory_groups
export consumer_groups
export prey_groups
export default_plankton_dynamics
export default_community
export default_biogeochem_dynamics

"""Return the fixed group set for a factory.

The returned tuple defines the canonical group order for group-level parameters.
"""
function factory_groups(::AbstractBGCFactory)
    throw(ArgumentError("No method `factory_groups(factory)` is defined for this factory."))
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

"""Default plankton community structure for a factory.

Returns a `NamedTuple` mapping group symbols to group specifications.

This is structural information only (group symbols, diameter specifications,
PFT specifications, etc.). Numeric parameter defaults are sourced from the
factory's default parameter generator (`Constructor.default_parameters`).
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

end # module FactoryInterface
