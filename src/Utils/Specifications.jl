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
export BiogeochemistrySpecification
export ModelSpecification

# -----------------------------------------------------------------------------
# Flexible specification containers
# -----------------------------------------------------------------------------

"""
    PFTSpecification(; kwargs...)

Container for plankton functional-type (PFT) specification parameters.

Stores arbitrary fields in `pft.data` (typically a `NamedTuple`).
"""
struct PFTSpecification
    data::Any
end

PFTSpecification(; kwargs...) = PFTSpecification((; kwargs...))

@inline pft_has(pft::PFTSpecification, key::Symbol) = hasproperty(pft.data, key)

@inline function pft_get(pft::PFTSpecification, key::Symbol, default=nothing)
    return pft_has(pft, key) ? getproperty(pft.data, key) : default
end

"""
    BiogeochemistrySpecification(; kwargs...)

Container for non-plankton biogeochemical specification parameters.

Stores arbitrary fields in `spec.data` (typically a `NamedTuple`).
"""
struct BiogeochemistrySpecification
    data::Any
end

BiogeochemistrySpecification(; kwargs...) = BiogeochemistrySpecification((; kwargs...))

"""
    ModelSpecification(data::NamedTuple)

Runtime parameter/specification container used by the generated biogeochemistry types.

This struct is `Adapt.jl`-compatible, so models can be adapted to GPU arrays.
"""
struct ModelSpecification{NT<:NamedTuple}
    data::NT
end

Adapt.@adapt_structure ModelSpecification

end # module Specifications
