module Models

using Agate.Utils: AbstractBGCFactory

# -----------------------------------------------------------------------------
# Factory default interface
# -----------------------------------------------------------------------------

"""Default plankton dynamics for a factory.

Returns a `NamedTuple` mapping a group prefix (e.g., `:P`, `:Z`) to a tracer
Dynamics builder function.
"""
function default_plankton_dynamics(::AbstractBGCFactory)
    throw(ArgumentError("No method `default_plankton_dynamics(factory)` is defined for this factory."))
end

"""Default plankton community structure for a factory.

Returns a `NamedTuple` mapping group prefix symbols to group specifications.

This is **structural** information only (group symbols, diameter specifications,
PFT specifications, etc.). All numeric parameter defaults live exclusively in
`Parameters.parameter_registry(factory)`.
"""
function default_community(::AbstractBGCFactory)
    throw(ArgumentError("No method `default_community(factory)` is defined for this factory."))
end

"""Default non-plankton tracer dynamics for a factory.

Returns a `NamedTuple` mapping tracer symbols (e.g., `:N`, `:DIC`) to tracer
Dynamics builder functions.
"""
function default_biogeochem_dynamics(::AbstractBGCFactory)
    throw(ArgumentError("No method `default_biogeochem_dynamics(factory)` is defined for this factory."))
end

# -----------------------------------------------------------------------------
# Submodules and public constructor
# -----------------------------------------------------------------------------

include("NiPiZD/NiPiZD.jl")
include("DARWIN/DARWIN.jl")

using .NiPiZD: NiPiZDFactory
using .DARWIN: DarwinFactory
export NiPiZDFactory
export DarwinFactory

end # module
