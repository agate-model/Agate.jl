"""Abstract supertype for canonical diameter specifications."""
abstract type AbstractDiameterSpecification end

"""A canonical diameter specification defined by explicit class diameters."""
struct DiameterListSpecification{V<:AbstractVector} <: AbstractDiameterSpecification
    diameters::V
end

"""A canonical diameter specification defined by a class count, bounds, and spacing rule."""
struct DiameterRangeSpecification{I<:Integer,T1,T2} <: AbstractDiameterSpecification
    n::I
    min_diameter::T1
    max_diameter::T2
    splitting::Symbol
end

function _validate_diameter_values(values, path::AbstractString)
    isempty(values) && throw(ArgumentError("$path must define at least one size class"))
    for (i, value) in pairs(values)
        value isa Real && !(value isa Bool) ||
            throw(ArgumentError("$path[$i] must be a real number"))
        isfinite(value) || throw(ArgumentError("$path[$i] must be finite; got $value"))
        value > zero(value) || throw(ArgumentError("$path[$i] must be positive; got $value"))
    end
    return nothing
end

"""Canonicalize and validate one public diameter input.

Accepted forms are an explicit vector, a `DiameterListSpecification`, a
`DiameterRangeSpecification`, or `(n, min_esd, max_esd, splitting)` as a NamedTuple.
"""
function normalize_diameters(diameters::AbstractVector; path::AbstractString="diameters")
    _validate_diameter_values(diameters, path)
    return (; n=length(diameters), specification=DiameterListSpecification(diameters))
end

function normalize_diameters(spec::NamedTuple; path::AbstractString="diameters")
    required = (:n, :min_esd, :max_esd, :splitting)
    all(hasproperty(spec, field) for field in required) || throw(
        ArgumentError("$path must define `n`, `min_esd`, `max_esd`, and `splitting`"),
    )
    return _normalize_diameter_range(
        spec.n, spec.min_esd, spec.max_esd, spec.splitting; path
    )
end

function normalize_diameters(
    spec::DiameterListSpecification; path::AbstractString="diameters"
)
    _validate_diameter_values(spec.diameters, path)
    return (; n=length(spec.diameters), specification=spec)
end

function normalize_diameters(
    spec::DiameterRangeSpecification; path::AbstractString="diameters"
)
    return _normalize_diameter_range(
        spec.n, spec.min_diameter, spec.max_diameter, spec.splitting; path
    )
end

normalize_diameters(spec; path::AbstractString="diameters") =
    throw(ArgumentError("$path has an unsupported diameter specification $(typeof(spec))"))

function _normalize_diameter_range(n, min_diameter, max_diameter, splitting; path)
    n isa Integer && !(n isa Bool) && n > 0 ||
        throw(ArgumentError("$path.n must be a positive integer; got $n"))
    min_diameter isa Real && !(min_diameter isa Bool) ||
        throw(ArgumentError("$path.min_esd must be a real number"))
    max_diameter isa Real && !(max_diameter isa Bool) ||
        throw(ArgumentError("$path.max_esd must be a real number"))
    isfinite(min_diameter) || throw(
        ArgumentError("$path.min_esd must be finite; got $min_diameter")
    )
    isfinite(max_diameter) || throw(
        ArgumentError("$path.max_esd must be finite; got $max_diameter")
    )
    min_diameter > zero(min_diameter) || throw(
        ArgumentError("$path.min_esd must be positive; got $min_diameter")
    )
    max_diameter > zero(max_diameter) || throw(
        ArgumentError("$path.max_esd must be positive; got $max_diameter")
    )
    max_diameter >= min_diameter || throw(
        ArgumentError("$path requires min_esd <= max_esd; got $min_diameter > $max_diameter")
    )
    splitting in (:log_splitting, :linear_splitting) || throw(
        ArgumentError(
            "$path.splitting must be :log_splitting or :linear_splitting; got $(repr(splitting))",
        ),
    )
    return (;
        n=Int(n),
        specification=DiameterRangeSpecification(
            Int(n), min_diameter, max_diameter, splitting
        ),
    )
end

function param_compute_diameters(
    ::Type{T}, spec::DiameterListSpecification
) where {T<:Real}
    return T[T(value) for value in spec.diameters]
end

function param_compute_diameters(
    ::Type{T}, spec::DiameterRangeSpecification
) where {T<:Real}
    n = Int(spec.n)
    min_d = T(spec.min_diameter)
    max_d = T(spec.max_diameter)
    n == 1 && return T[min_d]

    if spec.splitting === :log_splitting
        log_min = log(min_d)
        step = (log(max_d) - log_min) / T(n - 1)
        return T[exp(log_min + T(i - 1) * step) for i in 1:n]
    end

    step = (max_d - min_d) / T(n - 1)
    return T[min_d + T(i - 1) * step for i in 1:n]
end
