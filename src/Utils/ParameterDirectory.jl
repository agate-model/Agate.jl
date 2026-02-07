"""Parameter metadata for safer configuration.

Factories define `parameter_directory` entries that declare each required
parameter key, its shape (`:scalar`, `:vector`, or `:matrix`), and optional axes
for matrices.

The directory is used for validation and clearer error messages during
construction.
"""

export ParameterSpec
export DefaultProvider
export ParameterDefinition
export parameter_definitions
export ConstDefault
export NoDefault
export FillDefault
export DiameterIndexedVectorDefault
export parameter_directory
export parameter_spec

"""A small metadata record describing a model parameter."""
struct ParameterSpec
    name::Symbol
    shape::Symbol  # :scalar | :vector | :matrix
    axes::Union{Nothing,NTuple{2,Symbol}}
    doc::String
end

"""Convenience constructor for `ParameterSpec`."""
ParameterSpec(
    name::Symbol,
    shape::Symbol;
    axes::Union{Nothing,NTuple{2,Symbol}}=nothing,
    doc::AbstractString="",
) = ParameterSpec(name, shape, axes, String(doc))

# -----------------------------------------------------------------------------
# Single-source parameter definitions (spec + default)
# -----------------------------------------------------------------------------

"""Abstract supertype for constructor-time default providers.

Default providers are evaluated on the host during model construction.
They must produce concrete numeric values (scalars, vectors, matrices) that can
later be moved to a GPU architecture via `Adapt`.
"""
abstract type DefaultProvider end

"""Parameter definition that pairs a `ParameterSpec` with a default provider."""
struct ParameterDefinition{D<:DefaultProvider}
    spec::ParameterSpec
    default::D
end

"""A scalar default that converts a literal value to `FT`."""
struct ConstDefault{T} <: DefaultProvider
    value::T
end

"""Indicates that a parameter has no direct default value.

This is useful for parameters that are derived later (for example interaction
matrices regenerated from trait vectors).
"""
struct NoDefault <: DefaultProvider end

"""A uniform fill default for vectors or matrices."""
struct FillDefault{T} <: DefaultProvider
    value::T
end

"""Default provider for vectors defined over a subset of diameter-indexed classes.

The provider fills a full-length vector (length `community_context.n_total`) with
`default`, then overwrites the indices stored in `indices_field` (a field of
`CommunityContext`, e.g. `:producer_param_indices`) using
`resolve_diameter_indexed_vector`.
"""
struct DiameterIndexedVectorDefault{V,T} <: DefaultProvider
    value::V
    indices_field::Symbol
    default::T
end

DiameterIndexedVectorDefault(value, indices_field::Symbol; default) =
    DiameterIndexedVectorDefault(value, indices_field, default)

"""Return a tuple of `ParameterDefinition` entries for `factory`.

Factories should define one entry for every parameter key required by their
compiled equations.
"""
parameter_definitions(::AbstractBGCFactory) = ()

"""Return a tuple of `ParameterSpec` entries for `factory`.

By default the directory is derived from `parameter_definitions(factory)`.
Factories may still overload `parameter_directory` directly if needed.
"""
parameter_directory(factory::AbstractBGCFactory) = map(d -> d.spec, parameter_definitions(factory))

"""Return the `ParameterSpec` for `key`, or `nothing` if absent."""
function parameter_spec(factory::AbstractBGCFactory, key::Symbol)
    for spec in parameter_directory(factory)
        spec.name === key && return spec
    end
    return nothing
end
