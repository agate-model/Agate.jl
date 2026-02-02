"""Parameter metadata for safer configuration.

Factories define `parameter_directory` entries that declare each required
parameter key, its shape (`:scalar`, `:vector`, or `:matrix`), and optional axes
for matrices.

The directory is used for validation and clearer error messages during
construction.
"""

export ParameterSpec
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
