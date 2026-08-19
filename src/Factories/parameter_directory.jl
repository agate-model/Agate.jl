export ParameterSpec
export ParameterProvision
export DefaultProvider
export ParameterDefinition
export parameter_definitions
export ConstDefault
export NoDefault
export FillDefault
export DiameterIndexedVectorDefault
export DiameterIndexedMaterialization
export parameter_directory
export parameter_spec

"""Materialize a parameter law over selected diameter-indexed classes.

`role` names a declared parameter-applicability role in the parsed community.
When it is `nothing`, the law applies to every plankton class. `fill_value` is
assigned outside the selected role.
"""
struct DiameterIndexedMaterialization{T}
    role::Union{Nothing,Symbol}
    fill_value::T
end

function DiameterIndexedMaterialization(role::Union{Nothing,Symbol}=nothing; fill_value)
    return DiameterIndexedMaterialization(role, fill_value)
end

"""Declare that a model parameter supplies one semantic formulation requirement.

`process` is the stable named process identity, `path` locates an optional nested
sub-formulation, and `slot` is the formulation-local semantic parameter name.
`qualifier` distinguishes repeated instances of the same slot, such as Monod `K`
for different resources.
"""
struct ParameterProvision{P<:Tuple,Q<:NamedTuple}
    process::Symbol
    path::P
    formulation::Symbol
    slot::Symbol
    qualifier::Q
end

function ParameterProvision(
    process::Symbol,
    path::Tuple,
    formulation::Symbol,
    slot::Symbol;
    qualifier::NamedTuple=NamedTuple(),
)
    all(item -> item isa Symbol, path) || throw(
        ArgumentError("parameter provision path must contain only Symbols"),
    )
    return ParameterProvision(process, path, formulation, slot, qualifier)
end

function _parameter_provisions(provides)
    isnothing(provides) && return ()
    provides isa ParameterProvision && return (provides,)
    provides isa Tuple || throw(
        ArgumentError("`provides` must be a ParameterProvision or tuple of provisions"),
    )
    all(provision -> provision isa ParameterProvision, provides) || throw(
        ArgumentError("`provides` must contain only ParameterProvision values"),
    )
    return provides
end

"""Describe a configurable model parameter.

Fields
------
- `name`: parameter key.
- `shape`: one of `:scalar`, `:vector`, or `:matrix`.
- `axes`: optional legacy runtime vector axis name or matrix-axis names.
- `materialization`: optional constructor-time parameter-law materialization semantics.
- `provides`: semantic process/formulation requirements supplied by this parameter.
"""
struct ParameterSpec{P<:Tuple}
    name::Symbol
    shape::Symbol
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}
    materialization::Union{Nothing,DiameterIndexedMaterialization}
    provides::P
end

ParameterSpec(
    name::Symbol,
    shape::Symbol,
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}},
    materialization::Union{Nothing,DiameterIndexedMaterialization},
) = ParameterSpec(name, shape, axes, materialization, ())

"""Convenience constructor for `ParameterSpec`."""
function ParameterSpec(
    name::Symbol,
    shape::Symbol;
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}=nothing,
    materialization::Union{Nothing,DiameterIndexedMaterialization}=nothing,
    provides=(),
)
    provisions = _parameter_provisions(provides)
    return ParameterSpec(name, shape, axes, materialization, provisions)
end

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

"""A scalar default that converts a literal value to the construction scalar type."""
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
`default`, then overwrites the classes selected by the declared parameter role
using `resolve_diameter_indexed_vector`.
"""
struct DiameterIndexedVectorDefault{V,T} <: DefaultProvider
    value::V
    role::Symbol
    default::T
end

function DiameterIndexedVectorDefault(value, role::Symbol; default)
    return DiameterIndexedVectorDefault(value, role, default)
end

"""Return a tuple of `ParameterDefinition` entries for `factory`.

Factories should define one entry for every parameter key required by their
compiled equations.
"""
parameter_definitions(::AbstractBGCFactory) = ()

"""Return a tuple of `ParameterSpec` entries for `factory`.

By default the directory is derived from `parameter_definitions(factory)`.
Factories may still overload `parameter_directory` directly if needed.
"""
function parameter_directory(factory::AbstractBGCFactory)
    map(d -> d.spec, parameter_definitions(factory))
end

"""Return the `ParameterSpec` for `key`, or `nothing` if absent."""
function parameter_spec(factory::AbstractBGCFactory, key::Symbol)
    for spec in parameter_directory(factory)
        spec.name === key && return spec
    end
    return nothing
end
