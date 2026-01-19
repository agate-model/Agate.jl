"""Parameter metadata for safer configuration.

Agate's public constructors accept `NamedTuple` overrides for parameters and
interaction-related matrices. To keep this user-facing configuration explicit
(and to avoid brittle naming conventions), each factory can provide a
`parameter_directory` describing the expected *shape* and *kind* of each
parameter key.

The directory is small by design: it is intended for validation, better error
messages, and lightweight introspection during construction.
"""

export ParameterSpec
export parameter_directory
export parameter_spec

"""A small metadata record describing a model parameter."""
struct ParameterSpec
    name::Symbol
    shape::Symbol  # :scalar | :vector | :matrix
    kind::Symbol   # :real | :bool (introspection / error messages)
    doc::String
end

"""Convenience constructor for `ParameterSpec`."""
ParameterSpec(name::Symbol, shape::Symbol; kind::Symbol=:real, doc::AbstractString="") =
    ParameterSpec(name, shape, kind, String(doc))

"""Return a tuple of `ParameterSpec` entries for `factory`.

Factories should provide one entry for every parameter key required by their
compiled equations.
"""
parameter_directory(::AbstractBGCFactory) = ()

"""Return the `ParameterSpec` for `key`, or `nothing` if absent."""
function parameter_spec(factory::AbstractBGCFactory, key::Symbol)
    for spec in parameter_directory(factory)
        spec.name === key && return spec
    end
    return nothing
end
