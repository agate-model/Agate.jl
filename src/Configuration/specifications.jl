# """Agate.Configuration
# 
# Specification/container types used by Agate factories and constructors.
# 
# These are lightweight wrappers around `NamedTuple`s that support:
# 
# - ergonomic keyword-based overrides,
# - consistent casting to a target float type `FT` for CPU/GPU execution,
# - `Adapt.jl` compatibility for runtime parameter structs.
# 
# These types live under `Agate.Configuration` and are re-exported from
# `Agate.Construction` as part of the public factory/constructor API.
# """
using Adapt

export PFTSpecification, pft_get, pft_has

# -----------------------------------------------------------------------------
# Flexible specification containers
# -----------------------------------------------------------------------------

"""
    PFTSpecification(; kwargs...)

Container for plankton functional-type (PFT) specifications.

This is a lightweight wrapper around a `NamedTuple` used to attach per-PFT traits
and metadata.
"""
struct PFTSpecification
    data::Any
end

PFTSpecification(; kwargs...) = PFTSpecification((; kwargs...))

@inline pft_has(pft::PFTSpecification, key::Symbol) = hasproperty(pft.data, key)

@inline function pft_get(pft::PFTSpecification, key::Symbol, default=nothing)
    return pft_has(pft, key) ? getproperty(pft.data, key) : default
end

# Adapt support for `NamedTuple` payloads
#
# Many Agate runtime bundles store parameters in `NamedTuple`s for cheap field
# access and compile-time key sets. We rely on `Adapt.jl` to move these payloads
# to the architecture's preferred array type.
#
# `Adapt.jl` defines recursion for `NamedTuple` in modern releases. Defining it
# unconditionally overwrites upstream methods and triggers warnings during docs
# builds / precompilation. We only define a fallback rule when needed.
if !hasmethod(Adapt.adapt_structure, Tuple{Any,NamedTuple})
    @inline function Adapt.adapt_structure(to, nt::NamedTuple{names}) where {names}
        return NamedTuple{names}(map(x -> Adapt.adapt(to, x), values(nt)))
    end
end
