export ParameterSpec
export ParameterProvision
export DefaultProvider
export ParameterDefinition
export parameter_definitions
export ConstDefault
export DerivedDefault
export NoDefault
export FillDefault
export DiameterIndexedVectorDefault
export DiameterIndexedMaterialization
export derive_default
export parameter_directory
export parameter_spec

"""Materialize a diameter-indexed parameter law during construction.

Process participation determines which concrete classes the parameter applies to;
`fill_value` is assigned outside that applicability set.
"""
struct DiameterIndexedMaterialization{T}
    fill_value::T
end

DiameterIndexedMaterialization(; fill_value) = DiameterIndexedMaterialization(fill_value)

"""Declare that a model parameter supplies one semantic process parameter slot.

`process` is the stable named process identity and `slot` is the scientific parameter
name. `qualifier` narrows repeated slots, such as Monod `K` for a specific resource.
`path` is an optional disambiguator for models that contain more than one matching
slot within the same process. Formulation and resolved path are otherwise derived
from the normalized process definition.
"""
struct ParameterProvision{P<:Union{Nothing,Tuple},Q<:NamedTuple}
    process::Symbol
    slot::Symbol
    path::P
    qualifier::Q
end

function ParameterProvision(
    process::Symbol,
    slot::Symbol;
    path::Union{Nothing,Tuple}=nothing,
    qualifier::NamedTuple=NamedTuple(),
)
    isnothing(path) || all(item -> item isa Symbol, path) || throw(
        ArgumentError("parameter provision path must contain only Symbols"),
    )
    return ParameterProvision(process, slot, path, qualifier)
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
- `axes`: optional runtime vector axis name or matrix-axis names.
- `runtime_path`: concrete runtime storage path below `bgc.parameters`.
- `materialization`: optional constructor-time parameter-law materialization semantics.
- `provides`: semantic process parameter slots supplied by this parameter.
"""
struct ParameterSpec{R<:Tuple,P<:Tuple}
    name::Symbol
    shape::Symbol
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}
    runtime_path::R
    materialization::Union{Nothing,DiameterIndexedMaterialization}
    provides::P
end

"""Convenience constructor for `ParameterSpec`."""
function ParameterSpec(
    name::Symbol,
    shape::Symbol;
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}=nothing,
    runtime_path::Tuple=(name,),
    materialization::Union{Nothing,DiameterIndexedMaterialization}=nothing,
    provides=(),
)
    isempty(runtime_path) && throw(ArgumentError("runtime_path cannot be empty"))
    all(item -> item isa Symbol, runtime_path) || throw(
        ArgumentError("runtime_path must contain only Symbols"),
    )
    provisions = _parameter_provisions(provides)
    return ParameterSpec(name, shape, axes, runtime_path, materialization, provisions)
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

"""Define a model parameter and its constructor-time default.

`shape` is `:scalar`, `:vector`, or `:matrix`; vector and matrix definitions may
declare runtime `axes`. `provides` links the parameter to one or more scientific
process parameter slots.
"""
function ParameterDefinition(
    name::Symbol,
    default::D;
    shape::Symbol=:scalar,
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}=nothing,
    runtime_path::Tuple=(name,),
    materialization::Union{Nothing,DiameterIndexedMaterialization}=nothing,
    provides=(),
) where {D<:DefaultProvider}
    spec = ParameterSpec(name, shape; axes, runtime_path, materialization, provides)
    return ParameterDefinition(spec, default)
end

"""A scalar default that converts a literal value to the construction scalar type."""
struct ConstDefault{T} <: DefaultProvider
    value::T
end

"""Derive a parameter default from other resolved parameters during construction.

`deriver` is a setup-time strategy object implementing [`derive_default`](@ref).
`deps` names the parameter definitions that must be resolved before this default and
whose explicit overrides trigger recomputation.
"""
struct DerivedDefault{D,P<:Tuple} <: DefaultProvider
    deriver::D
    deps::P
end

function DerivedDefault(deriver; deps=())
    deps isa Tuple || throw(ArgumentError("deps must be a tuple of Symbols."))
    all(dep -> dep isa Symbol, deps) || throw(
        ArgumentError("deps must be a tuple of Symbols.")
    )
    length(Set(deps)) == length(deps) || throw(
        ArgumentError("deps must not contain duplicate parameter keys.")
    )
    return DerivedDefault{typeof(deriver),typeof(deps)}(deriver, deps)
end

"""Compute a value for a [`DerivedDefault`](@ref) provider.

Concrete derivers receive the owning model source, construction context, and all
parameter values resolved so far. The owner is a registered model family for named
models and the authored `ModelDefinition` for direct construction. Derivation runs on
the host during model construction.
"""
function derive_default(deriver, owner, context, params::NamedTuple)
    throw(ArgumentError("derive_default is not implemented for $(typeof(deriver))."))
end

"""Indicates that a parameter has no default value.

The parameter must be supplied by an override before construction can complete.
"""
struct NoDefault <: DefaultProvider end

"""A uniform fill default for vectors or matrices."""
struct FillDefault{T} <: DefaultProvider
    value::T
end

"""Default provider for a diameter-indexed vector parameter.

The provider fills the complete runtime vector with `default`, then materializes
`value` over the classes selected from process-derived parameter applicability.
"""
struct DiameterIndexedVectorDefault{V,T} <: DefaultProvider
    value::V
    default::T
end

DiameterIndexedVectorDefault(value; default) = DiameterIndexedVectorDefault(value, default)

"""Return the parameter definitions owned by a model family or direct definition source."""
parameter_definitions(source) = ()

"""Return the parameter specifications derived from `parameter_definitions(source)`."""
parameter_directory(source) = map(d -> d.spec, parameter_definitions(source))

"""Return the `ParameterSpec` for `key`, or `nothing` if absent."""
function parameter_spec(source, key::Symbol)
    for spec in parameter_directory(source)
        spec.name === key && return spec
    end
    return nothing
end
