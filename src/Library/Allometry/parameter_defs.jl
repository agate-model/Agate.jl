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

    # By construction we keep coefficients and diameter the same floating type (FT),
    # so this call never mixes Float32/Float64 (important for GPU use).
    return allometric_scaling_power(a, b, diameter)
end

"""Resolve a parameter definition at a given diameter.

This is the one function the constructor uses when building runtime parameter vectors.
"""
@inline resolve_param(::Type{FT}, x, diameter) where {FT<:AbstractFloat} = FT(x)

@inline resolve_param(::Type{FT}, x::Bool, diameter) where {FT<:AbstractFloat} = x

@inline resolve_param(::Type{FT}, p::ConstantParam, diameter) where {FT<:AbstractFloat} =
    FT(p.value)

@inline function resolve_param(
    ::Type{FT}, p::AllometricParam, diameter
) where {FT<:AbstractFloat}
    # Coefficients often come from literal numbers (Float64). Convert them to FT so we
    # never mix Float32/Float64 in the underlying allometric calls.
    coeffs = map(v -> v isa Number ? FT(v) : v, p.coeffs)
    return FT(p.model(coeffs, FT(diameter)))
end

"""Resolve a parameter definition across a set of diameters.

Returns a vector `v` where `v[i] = resolve_param(FT, value, diameters[i])`.

This helper is typically used when building constructor-time default parameter
vectors from scalar or allometric definitions.
"""
function resolve_diameter_vector(
    ::Type{FT}, diameters::AbstractVector, value
) where {FT<:AbstractFloat}
    n = length(diameters)
    out = Vector{FT}(undef, n)
    @inbounds for i in 1:n
        out[i] = resolve_param(FT, value, diameters[i])
    end
    return out
end

"""Resolve a parameter definition over a subset of diameters.

Returns a vector of length `length(diameters)` filled with `default`, then
overwrites entries indexed by `indices` with `resolve_param(FT, value, diameters[i])`.

This is convenient for parameters that apply only to certain plankton roles
(e.g. producer-only growth rates or consumer-only predation rates).
"""
function resolve_diameter_indexed_vector(
    ::Type{FT},
    diameters::AbstractVector,
    indices::AbstractVector{<:Integer},
    value;
    default::FT,
) where {FT<:AbstractFloat}
    out = fill(default, length(diameters))
    @inbounds for i in indices
        out[i] = resolve_param(FT, value, diameters[i])
    end
    return out
end

"""Cast numeric entries inside a parameter definition to `FT`."""
@inline function cast_paramdef(::Type{FT}, p::ConstantParam) where {FT<:AbstractFloat}
    return ConstantParam(FT(p.value))
end

@inline function cast_paramdef(::Type{FT}, p::AllometricParam) where {FT<:AbstractFloat}
    coeffs = p.coeffs
    if coeffs isa NamedTuple
        coeffs = map(v -> v isa Number ? FT(v) : v, coeffs)
    end
    return AllometricParam(p.model, coeffs)
end
