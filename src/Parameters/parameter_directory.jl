export DefaultProvider
export Parameter
export parameter_definitions
export DerivedDefault
export NoDefault
export DiameterIndexedVectorDefault
export derive_default
export parameter_directory
export parameter_spec

"""Describe one configurable model parameter independent of its model-level name.

The parameter name is owned by the key in the model's `parameters` NamedTuple. `axes`
optionally selects explicit runtime storage axes. With `axes=nothing`, dimensionality is
resolved structurally from the bound scientific slot; dependency-only parameters are scalar.
"""
struct ParameterSpec
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}
end

"""Convenience constructor for `ParameterSpec`."""
ParameterSpec(; axes::Union{Nothing,Symbol,NTuple{2,Symbol}}=nothing) = ParameterSpec(axes)

"""Abstract supertype for constructor-time default providers.

Default providers are evaluated on the host during model construction.
They must produce concrete numeric values (scalars, vectors, matrices) that can
later be moved to a GPU architecture via `Adapt`.
"""
abstract type DefaultProvider end

"""A keyed model parameter with setup-time metadata and a default provider.

The stable model parameter name is the key in the enclosing `parameters` NamedTuple,
for example `parameters=(maximum_growth_rate=Parameter(...),)`. `Parameter` therefore
does not duplicate its name internally.
"""
struct Parameter{D<:DefaultProvider}
    spec::ParameterSpec
    default::D
end

"""Define one model parameter value/default independent of its model-level key.

Explicit `axes` select vector (`Symbol`) or matrix (`NTuple{2,Symbol}`) runtime storage.
With `axes=nothing`, dimensionality comes from the bound scientific slot; a parameter used
only as a `DerivedDefault` dependency is scalar.
"""
function Parameter(
    default::D;
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}=nothing,
) where {D<:DefaultProvider}
    return Parameter(ParameterSpec(; axes), default)
end

"""Define a parameter with a literal constant default.

Literal defaults are wrapped in `ConstantDefault` automatically; explicit default-provider
objects remain available for derived, missing, or structure-dependent defaults.
"""
function Parameter(
    default;
    axes::Union{Nothing,Symbol,NTuple{2,Symbol}}=nothing,
)
    return Parameter(ConstantDefault(default); axes)
end

"""A literal parameter default.

Scalar parameters use the literal value directly after conversion to the construction
scalar type. Vector and matrix parameters fill their resolved storage extent with the
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
`value` over process-derived applicability. The same `default` value is used outside
applicability when a diameter-indexed parameter-law override is supplied. Parameters without process applicability
materialize over the full living-class axis.
"""
struct DiameterIndexedVectorDefault{V,T} <: DefaultProvider
    value::V
    default::T
end

DiameterIndexedVectorDefault(value; default) = DiameterIndexedVectorDefault(value, default)

"""Return the keyed parameter definitions owned by a model family or direct source."""
parameter_definitions(source) = (;)

"""Return authored parameter specifications keyed by stable model parameter name."""
function parameter_directory(source)
    definitions = parameter_definitions(source)
    definitions isa NamedTuple || throw(
        ArgumentError("parameter_definitions(::$(typeof(source))) must return a NamedTuple"),
    )
    return NamedTuple{keys(definitions)}(Tuple(definition.spec for definition in values(definitions)))
end

"""Return the `ParameterSpec` for `key`, or `nothing` if absent."""
function parameter_spec(source, key::Symbol)
    directory = parameter_directory(source)
    return hasproperty(directory, key) ? getproperty(directory, key) : nothing
end
