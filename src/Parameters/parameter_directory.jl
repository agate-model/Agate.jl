export ParameterSpec
export ParameterProvision
export DefaultProvider
export ParameterDefinition
export parameter_definitions
export ConstantDefault
export DerivedDefault
export NoDefault
export DiameterIndexedVectorDefault
export DiameterIndexedMaterialization
export derive_default
export parameter_directory
export parameter_spec

"""Materialize a diameter-indexed parameter law during construction.

Process participation determines which concrete classes a provision-bound parameter applies
to; provision-less diameter-indexed parameters use the full living-class axis. `fill_value`
is assigned outside the selected applicability set.
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
- `shape`: resolved storage shape (`:scalar`, `:vector`, or `:matrix`); may be `nothing` before process-bound shape inference.
- `axes`: optional explicit runtime storage axis/axes. With `axes=nothing`,
  vector/matrix storage is process-local and follows the resolved provision applicability.
  `axes=:plankton` selects the full living-class axis; matrix axes such as
  `(:consumer, :prey)` select the corresponding global role axes.
- `materialization`: optional constructor-time parameter-law materialization semantics.
- `provides`: semantic process parameter slots supplied by this parameter.
"""
struct ParameterSpec{P<:Tuple}
    name::Symbol
    shape::Union{Nothing,Symbol}
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}
    materialization::Union{Nothing,DiameterIndexedMaterialization}
    provides::P
end

"""Convenience constructor for `ParameterSpec`."""
function ParameterSpec(
    name::Symbol,
    shape::Union{Nothing,Symbol}=nothing;
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}=nothing,
    materialization::Union{Nothing,DiameterIndexedMaterialization}=nothing,
    provides=(),
)
    provisions = _parameter_provisions(provides)
    declared_shape = if isnothing(shape)
        axes isa Symbol ? :vector : axes isa Tuple ? :matrix : nothing
    else
        shape
    end
    return ParameterSpec(name, declared_shape, axes, materialization, provisions)
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

`shape` may be `:scalar`, `:vector`, or `:matrix`. For provision-bound parameters it
may be omitted and is inferred from the resolved process requirement. Explicit runtime
`axes` also imply vector or matrix shape. Parameters with neither provisions nor axes
must declare shape explicitly. `provides` links the parameter to one or more scientific
process parameter slots.
"""
function ParameterDefinition(
    name::Symbol,
    default::D;
    shape::Union{Nothing,Symbol}=nothing,
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}=nothing,
    materialization::Union{Nothing,DiameterIndexedMaterialization}=nothing,
    provides=(),
) where {D<:DefaultProvider}
    spec = ParameterSpec(name, shape; axes, materialization, provides)
    return ParameterDefinition(spec, default)
end

"""A literal parameter default.

Scalar parameters use the literal value directly after conversion to the construction
scalar type. Vector and matrix parameters fill their resolved storage shape with the
literal value.
"""
struct ConstantDefault{T} <: DefaultProvider
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

"""Default provider for a diameter-indexed vector parameter.

The provider fills the complete runtime vector with `default`, then materializes
`value` over process-derived applicability. Provision-less parameters materialize over
the full living-class axis.
"""
struct DiameterIndexedVectorDefault{V,T} <: DefaultProvider
    value::V
    default::T
end

DiameterIndexedVectorDefault(value; default) = DiameterIndexedVectorDefault(value, default)

"""Return the parameter definitions owned by a model family or direct definition source."""
parameter_definitions(source) = ()

"""Return authored parameter specifications derived from `parameter_definitions(source)`.

Provision-bound specifications may have `shape === nothing` until `normalize_model` resolves
their process requirements.
"""
parameter_directory(source) = map(d -> d.spec, parameter_definitions(source))

"""Return the `ParameterSpec` for `key`, or `nothing` if absent."""
function parameter_spec(source, key::Symbol)
    for spec in parameter_directory(source)
        spec.name === key && return spec
    end
    return nothing
end
