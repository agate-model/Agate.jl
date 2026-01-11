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

export PFTSpecification, pft_get, pft_has, cast_pft
export BiogeochemistrySpecification, cast_spec
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

# -----------------------------------------------------------------------------
# Casting utilities
# -----------------------------------------------------------------------------

@inline _cast_number(::Type{FT}, x) where {FT<:AbstractFloat} = FT(x)

function _cast_container(::Type{FT}, x) where {FT<:AbstractFloat}
    if x isa Bool
        return x
    elseif x isa Number
        return _cast_number(FT, x)
    elseif x isa AbstractArray
        return x isa AbstractArray{Bool} ? x : FT.(x)
    elseif x isa NamedTuple
        return map(v -> _cast_container(FT, v), x)
    else
        return x
    end
end

"""Cast numeric entries in `pft` to `FT` (recursively for arrays and NamedTuples)."""
function cast_pft(::Type{FT}, pft::PFTSpecification) where {FT<:AbstractFloat}
    return PFTSpecification(_cast_container(FT, pft.data))
end

"""Cast numeric entries in `spec` to `FT` (recursively for arrays and NamedTuples)."""
function cast_spec(::Type{FT}, spec::BiogeochemistrySpecification) where {FT<:AbstractFloat}
    return BiogeochemistrySpecification(_cast_container(FT, spec.data))
end

end # module Specifications
