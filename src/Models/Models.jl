module Models

using Agate.Utils: AbstractBGCFactory

# -----------------------------------------------------------------------------
# Factory default interface
# -----------------------------------------------------------------------------

"""Default plankton dynamics for a factory.

Returns a `NamedTuple` mapping a group prefix (e.g., `:P`, `:Z`) to a tracer
dynamics builder function.
"""
function default_plankton_dynamics(::AbstractBGCFactory)
    throw(ArgumentError("No method `default_plankton_dynamics(factory)` is defined for this factory."))
end

"""Default plankton arguments for a factory.

Returns a `NamedTuple` mapping group prefix symbols to group specifications.
"""
function default_plankton_args(::AbstractBGCFactory, ::Type{FT}) where {FT<:AbstractFloat}
    throw(ArgumentError("No method `default_plankton_args(factory, FT)` is defined for this factory."))
end

"""Default non-plankton tracer dynamics for a factory.

Returns a `NamedTuple` mapping tracer symbols (e.g., `:N`, `:DIC`) to tracer
dynamics builder functions.
"""
function default_biogeochem_dynamics(::AbstractBGCFactory)
    throw(ArgumentError("No method `default_biogeochem_dynamics(factory)` is defined for this factory."))
end

"""Default biogeochemical specification for a factory."""
function default_biogeochem_args(::AbstractBGCFactory, ::Type{FT}) where {FT<:AbstractFloat}
    throw(ArgumentError("No method `default_biogeochem_args(factory, FT)` is defined for this factory."))
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
