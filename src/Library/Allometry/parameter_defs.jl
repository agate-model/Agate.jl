# -----------------------------------------------------------------------------
# Explicit parameter definitions
# -----------------------------------------------------------------------------

"""
    AbstractParamDef

Abstract supertype for constructor-time parameter definitions.

Concrete parameter definitions describe how scalar or size-dependent parameter
values are materialized during model construction.
"""
abstract type AbstractParamDef end

"""
    ConstantParam(value)

Parameter definition whose value is constant across all size classes.

!!! formulation
    ```math
    p(d_i) = p_0
    ```

    where ``p_0`` is `value` and ``d_i`` is any class diameter.
"""
struct ConstantParam{T} <: AbstractParamDef
    value::T
end

# NOTE: the default outer constructor `ConstantParam(x)` already exists and infers `T`.
# We avoid redefining it to prevent method overwrite warnings during precompilation.

"""
    AllometricParam(model, coeffs)
    AllometricParam(model; kwargs...)

Parameter definition evaluated from an allometric model and coefficient bundle.

!!! formulation
    ```math
    p(d_i) = model(coeffs, d_i)
    ```

    where `model` is a callable and `coeffs` is a named tuple of coefficients.
    The keyword constructor stores `kwargs` as the coefficient named tuple.
"""
struct AllometricParam{F,C} <: AbstractParamDef
    model::F
    coeffs::C
end

AllometricParam(model; kwargs...) = AllometricParam(model, (; kwargs...))

"""
    PowerLaw()

Callable allometric power-law model using spherical cell volume.

!!! formulation
    ```math
    p(d) = a V(d)^b,
    \\qquad
    V(d) = \\frac{4}{3}\\pi\\left(\\frac{d}{2}\\right)^3
    ```

    The expected coefficient names are `prefactor` for ``a`` and `exponent` for
    ``b``.
"""
struct PowerLaw end

"""
    PowerLaw()(coeffs, diameter)

Evaluate a `PowerLaw` allometric model.

!!! formulation
    ```math
    p(d) = a V(d)^b,
    \\qquad
    V(d) = \\frac{4}{3}\\pi\\left(\\frac{d}{2}\\right)^3
    ```

# Arguments
- `coeffs`: named tuple with `prefactor` and `exponent` entries.
- `diameter`: equivalent spherical diameter ``d``.
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

"""
    resolve_param(T, value, diameter)

Resolve a scalar, boolean, or explicit parameter definition at one diameter.

!!! formulation
    For scalar values, `resolve_param` returns `T(value)`. For
    `ConstantParam(p₀)`, it returns ``T(p_0)``. For
    `AllometricParam(model, coeffs)`, it returns

    ```math
    T\\left(model(coeffs_T, T(d))\\right)
    ```

    where numeric coefficients are converted to `T` before evaluation.

# Arguments
- `T`: target scalar type.
- `value`: scalar, boolean, `ConstantParam`, or `AllometricParam`.
- `diameter`: class diameter used by allometric definitions.
"""
@inline resolve_param(::Type{T}, x, diameter) where {T<:Real} = T(x)

@inline resolve_param(::Type{T}, x::Bool, diameter) where {T<:Real} = x

@inline resolve_param(::Type{T}, p::ConstantParam, diameter) where {T<:Real} = T(p.value)

@inline function resolve_param(::Type{T}, p::AllometricParam, diameter) where {T<:Real}
    # Coefficients often come from literal numbers (Float64). Convert them to the construction scalar type so we
    # never mix Float32/Float64 in the underlying allometric calls.
    coeffs = map(v -> v isa Number ? T(v) : v, p.coeffs)
    return T(p.model(coeffs, T(diameter)))
end

"""
    resolve_diameter_vector(T, diameters, value)

Resolve a scalar or parameter definition across all diameters.

!!! formulation
    ```math
    out_i = resolve\\_param(T, value, d_i)
    ```

    where ``d_i`` is `diameters[i]`.

# Arguments
- `T`: target scalar type.
- `diameters`: diameter vector.
- `value`: scalar or parameter definition to resolve.
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

"""
    resolve_diameter_indexed_vector(T, diameters, indices, value; default)

Resolve a scalar or parameter definition over a subset of diameters.

!!! formulation
    ```math
    out_i = \\begin{cases}
        resolve\\_param(T, value, d_i), & i \\in indices \\
        default, & otherwise
    \\end{cases}
    ```

    This is useful for parameters that apply only to selected model classes,
    such as producer-only growth rates or consumer-only grazing rates.

# Arguments
- `T`: target scalar type.
- `diameters`: diameter vector.
- `indices`: indices to resolve from `value`.
- `value`: scalar or parameter definition to resolve at selected indices.
- `default`: value assigned outside `indices`.
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

"""
    cast_paramdef(T, p)

Cast numeric entries inside a parameter definition to scalar type `T`.

!!! formulation
    `ConstantParam(p₀)` becomes `ConstantParam(T(p₀))`. For
    `AllometricParam(model, coeffs)`, numeric entries in `coeffs` are converted
    to `T` and non-numeric entries are preserved.
"""
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
