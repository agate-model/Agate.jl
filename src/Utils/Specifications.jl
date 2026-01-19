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

# Ensure `PFTSpecification` payloads adapt correctly to the target architecture.
#
# NOTE: Do *not* define `Adapt.adapt_structure` for `NamedTuple` here.
# `Adapt.jl` already provides a rule for `NamedTuple` and defining a second one
# triggers method overwrite warnings during documentation builds.
Adapt.@adapt_structure PFTSpecification

@inline pft_has(pft::PFTSpecification, key::Symbol) = hasproperty(pft.data, key)

@inline function pft_get(pft::PFTSpecification, key::Symbol, default=nothing)
    return pft_has(pft, key) ? getproperty(pft.data, key) : default
end

end # module Specifications
