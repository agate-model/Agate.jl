export DefaultProvider
export Parameter
export MetaParameter
export parameter_definitions
export DerivedDefault
export NoDefault
export DiameterIndexedVectorDefault
export derive_default

"""Abstract supertype for construction-time default providers.

Default providers are evaluated on the host during model construction.
They must produce concrete numeric values (scalars, vectors, matrices) that can
later be moved to a GPU architecture via `Adapt`.
"""
abstract type DefaultProvider end

"""A literal parameter default.

Scalar parameters use the literal value directly after conversion to the construction
scalar type. Vector and matrix parameters fill their resolved storage extent with the
literal value.
"""
struct ConstantDefault{T} <: DefaultProvider
    value::T
end

"""A runtime model parameter whose storage is determined by scientific process slots.

The stable parameter name is the key in the enclosing `parameters` NamedTuple. Scalar,
vector, and matrix storage are inferred from the slots that bind the parameter and their
realized process applicability; process slots therefore determine the storage axes.
"""
struct Parameter{D<:DefaultProvider}
    default::D
end

"""Define a runtime parameter with a literal constant default."""
Parameter(default) = Parameter(ConstantDefault(default))

"""A construction-only parameter used to derive one or more runtime parameters.

`MetaParameter` values are materialized during setup but are not stored in the runtime
biogeochemistry. Shaped meta-parameters use `axes=:plankton`; omitting `axes` defines a
scalar meta-parameter.
"""
struct MetaParameter{D<:DefaultProvider,A}
    default::D
    axes::A
end

function MetaParameter(
    default::D;
    axes=nothing,
) where {D<:DefaultProvider}
    (isnothing(axes) || axes === :plankton) || throw(
        ArgumentError("MetaParameter axes must be `nothing` or `:plankton`"),
    )
    return MetaParameter{D,typeof(axes)}(default, axes)
end

"""Define a construction-only parameter with a literal constant default."""
function MetaParameter(
    default;
    axes=nothing,
)
    return MetaParameter(ConstantDefault(default); axes)
end

"""Derive a parameter default from other resolved parameters during construction.

`deriver` is a setup-time strategy object implementing [`derive_default`](@ref).
`deps` names ordinary runtime parameters and/or construction-only [`MetaParameter`](@ref)
values available to the derivation. Derived defaults are evaluated once after direct
defaults and explicit overrides are materialized.
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

Concrete derivers receive the owning model source, construction context, and a
`NamedTuple` containing exactly the dependencies declared by `deps`. The owner is a
registered model family for named models and the authored `ModelDefinition` for direct
construction. Derivation runs on the host during model construction.
"""
function derive_default(deriver, owner, context, params::NamedTuple)
    throw(ArgumentError("derive_default is not implemented for $(typeof(deriver))."))
end

"""Indicates that a parameter has no default value.

The parameter must be supplied by an override before construction can complete.
"""
struct NoDefault <: DefaultProvider end

"""Default provider for a diameter-indexed vector parameter.

The provider materializes `value` over the parameter's realized storage labels. For
slot-bound [`Parameter`](@ref) values those labels come from process applicability; for a
[`MetaParameter`](@ref) they come from its explicit construction axis.
"""
struct DiameterIndexedVectorDefault{V,T} <: DefaultProvider
    value::V
    default::T
end

DiameterIndexedVectorDefault(value; default) = DiameterIndexedVectorDefault(value, default)

"""Return the keyed parameter definitions owned by a model family or direct source."""
parameter_definitions(source) = (;)
