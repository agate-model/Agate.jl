module Models

using ..Utils: AbstractBGCFactory

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

This is *structural* information only (group symbols, diameters, etc.).
All *parameter defaults* live exclusively in the model's parameter registry.
"""
function default_community(::AbstractBGCFactory, ::Type{FT}) where {FT<:AbstractFloat}
    throw(ArgumentError("No method `default_community(factory, FT)` is defined for this factory."))
end

"""Keyword front-end for `default_community(factory, FT)`.

The public constructors and tests pass the floating-point type as a keyword
argument (`FT=Float32`, `FT=Float64`). Internally we keep the canonical API as
`default_community(factory, ::Type{FT})` to preserve clean parametric
dispatch.
"""
function default_community(factory::AbstractBGCFactory; FT::Type{T}=Float64) where {T<:AbstractFloat}
    return default_community(factory, T)
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
