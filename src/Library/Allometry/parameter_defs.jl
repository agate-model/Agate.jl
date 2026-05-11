# -----------------------------------------------------------------------------
# Explicit parameter definitions
# -----------------------------------------------------------------------------

"""Abstract supertype for explicit parameter definitions stored in `PFTSpecification`."""
abstract type AbstractParamDef end

"""A parameter that is constant across sizes."""
struct ConstantParam{T} <: AbstractParamDef
    value::T
end

# NOTE: the default outer constructor `ConstantParam(x)` already exists and infers `T`.
# We avoid redefining it to prevent method overwrite warnings during precompilation.

"""A parameter that is computed from an explicit model/callable and coefficient bundle."""
struct AllometricParam{F,C} <: AbstractParamDef
    model::F
    coeffs::C
end

"""Construct an `AllometricParam` from a model and keyword coefficients."""
AllometricParam(model; kwargs...) = AllometricParam(model, (; kwargs...))

"""Common allometry model: power-law scaling with spherical volume."""
struct PowerLaw end

"""Evaluate a `PowerLaw` model.

Expected coefficient names:
- `prefactor`
- `exponent`

(Backward-compat coefficient aliases were intentionally removed.)
"""
@inline function (m::PowerLaw)(coeffs::NamedTuple, diameter)
    hasproperty(coeffs, :prefactor) ||
        throw(ArgumentError("PowerLaw requires coefficient `prefactor`"))
    hasproperty(coeffs, :exponent) ||
        throw(ArgumentError("PowerLaw requires coefficient `exponent`"))

    a = getproperty(coeffs, :prefactor)
    b = getproperty(coeffs, :exponent)

    # By construction we keep coefficients and diameter the same scalar type,
    # so this call never mixes Float32/Float64 (important for GPU use).
    return allometric_scaling_power(a, b, diameter)
end

"""Resolve a parameter definition at a given diameter.

This is the one function the constructor uses when building runtime parameter vectors.
"""
@inline resolve_param(::Type{T}, x, diameter) where {T<:Real} = T(x)

@inline resolve_param(::Type{T}, x::Bool, diameter) where {T<:Real} = x

@inline resolve_param(::Type{T}, p::ConstantParam, diameter) where {T<:Real} =
    T(p.value)

@inline function resolve_param(
    ::Type{T}, p::AllometricParam, diameter
) where {T<:Real}
    # Coefficients often come from literal numbers (Float64). Convert them to the construction scalar type so we
    # never mix Float32/Float64 in the underlying allometric calls.
    coeffs = map(v -> v isa Number ? T(v) : v, p.coeffs)
    return T(p.model(coeffs, T(diameter)))
end

"""Resolve a parameter definition across a set of diameters.

Returns a vector `v` where `v[i] = resolve_param(T, value, diameters[i])`.

This helper is typically used when building constructor-time default parameter
vectors from scalar or allometric definitions.
"""
function resolve_diameter_vector(
    ::Type{T}, diameters::AbstractVector, value
) where {T<:Real}
    n = length(diameters)
    out = Vector{T}(undef, n)
    @inbounds for i in 1:n
        out[i] = resolve_param(T, value, diameters[i])
    end
    return out
end

"""Resolve a parameter definition over a subset of diameters.

Returns a vector of length `length(diameters)` filled with `default`, then
overwrites entries indexed by `indices` with `resolve_param(T, value, diameters[i])`.

This is convenient for parameters that apply only to certain plankton roles
(e.g. producer-only growth rates or consumer-only predation rates).
"""
function resolve_diameter_indexed_vector(
    ::Type{T},
    diameters::AbstractVector,
    indices::AbstractVector{<:Integer},
    value;
    default::T,
) where {T<:Real}
    out = fill(default, length(diameters))
    @inbounds for i in indices
        out[i] = resolve_param(T, value, diameters[i])
    end
    return out
end

"""Cast numeric entries inside a parameter definition to `T`."""
@inline function cast_paramdef(::Type{T}, p::ConstantParam) where {T<:Real}
    return ConstantParam(T(p.value))
end

@inline function cast_paramdef(::Type{T}, p::AllometricParam) where {T<:Real}
    coeffs = p.coeffs
    if coeffs isa NamedTuple
        coeffs = map(v -> v isa Number ? T(v) : v, coeffs)
    end
    return AllometricParam(p.model, coeffs)
end
