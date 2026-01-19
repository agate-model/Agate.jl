"""Agate.Utils.Specifications

Specification/container types used by Agate factories and constructors.

These are lightweight wrappers around `NamedTuple`s that support:

- ergonomic keyword-based overrides,
- consistent casting to a target float type `FT` for CPU/GPU execution,
- `Adapt.jl` compatibility for runtime parameter structs.

These types live under `Agate.Utils.Specifications` and are re-exported from
`Agate.Constructor` as part of the public factory/constructor API.
"""
module Specifications

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
# `Adapt` does not always ship with a recursive rule for `NamedTuple` across all
# versions, so we define it here to ensure that `adapt(CuArray, params)` converts
# vectors/matrices inside the `NamedTuple`.
@inline function Adapt.adapt_structure(to, nt::NamedTuple{names}) where {names}
    return NamedTuple{names}(map(x -> Adapt.adapt(to, x), values(nt)))
end

end # module Specifications
