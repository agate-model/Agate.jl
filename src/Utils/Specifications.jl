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
export ModelSpecification

# -----------------------------------------------------------------------------
# Flexible specification containers
# -----------------------------------------------------------------------------

"""
    PFTSpecification(; kwargs...)

Container for plankton functional-type (PFT) specifications.

This is a lightweight wrapper around a `NamedTuple` used to attach per-PFT traits
and metadata (e.g. feeding traits used to build default interaction matrices).

The parameter registry system does **not** read `PFTSpecification` fields to override
registry parameters. Group-level parameters are configured via `update_registry`
(full replacement) and `patch_registry_groups` (explicit partial updates).
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
    ModelSpecification(data::NamedTuple)

Runtime parameter/specification container used by the generated biogeochemistry types.

This struct is `Adapt.jl`-compatible, so models can be adapted to GPU arrays.
"""
struct ModelSpecification{NT<:NamedTuple}
    data::NT
end

Adapt.@adapt_structure ModelSpecification

# -----------------------------------------------------------------------------
# Ergonomic access
# -----------------------------------------------------------------------------

"""Forward `getproperty` to the underlying `NamedTuple` payload.

This allows end users to write `params.foo` instead of `params.data.foo` while
still keeping `.data` available for introspection.
"""
@inline function Base.getproperty(ms::ModelSpecification, name::Symbol)
    name === :data && return getfield(ms, :data)
    return getproperty(getfield(ms, :data), name)
end

@inline Base.hasproperty(ms::ModelSpecification, name::Symbol) =
    (name === :data) || hasproperty(getfield(ms, :data), name)

@inline function Base.propertynames(ms::ModelSpecification, private::Bool=false)
    return (:data, propertynames(getfield(ms, :data), private)...)
end

end # module Specifications
