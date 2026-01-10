"""
Utilities to construct Oceananigans biogeochemistry types.

Agate constructs `AbstractContinuousFormBiogeochemistry` subtypes from:
- a runtime parameter struct stored as `bgc.parameters`, and
- tracer tendency expressions (`Expr`) keyed by tracer name.

The generated types are compatible with Oceananigans and OceanBioME.
Runtime parameter structs can be adapted between CPU and GPU with `Adapt.jl`.
"""
module Utils

using Adapt

using Agate.Library.Mortality
using Agate.Library.Nutrients
using Agate.Library.Photosynthesis
using Agate.Library.Predation
using Agate.Library.Remineralization

using OceanBioME
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Fields: ZeroField

import Oceananigans.Biogeochemistry:
    biogeochemical_drift_velocity,
    required_biogeochemical_auxiliary_fields,
    required_biogeochemical_tracers

export AbstractDiameterSpecification
export DiameterListSpecification
export DiameterRangeSpecification

export add_bgc_methods!, create_bgc_struct, define_tracer_functions, expression_check

export param_check_length
export param_check_square_matrix
export param_cast_matrix
export param_compute_diameters

#temp define here:
"""Abstract supertype for diameter specifications."""
abstract type AbstractDiameterSpecification end

"""A diameter specification defined by an explicit list of diameters."""
struct DiameterListSpecification{T,VT<:AbstractVector{T}} <: AbstractDiameterSpecification
    diameters::VT
end

"""A diameter specification defined by a range and a splitting method."""
struct DiameterRangeSpecification{T} <: AbstractDiameterSpecification
    min_diameter::T
    max_diameter::T
    splitting::Symbol
end

# -----------------------------------------------------------------------------
# Expression validation
# -----------------------------------------------------------------------------

"""Return all Symbols referenced by an expression tree."""
function parse_expression(f_expr)
    symbols = Symbol[]
    expressions = Expr[f_expr]

    for exp in expressions
        for arg in exp.args
            if arg isa Expr
                push!(expressions, arg)
            elseif arg isa Symbol
                push!(symbols, arg)
            end
        end
    end

    return symbols
end

"""
    expression_check(allowed_symbols, f_expr; module_name=Utils)

Validate that all Symbols referenced in `f_expr` are either:
- present in `allowed_symbols`, or
- defined in `module_name`.

Throws `UndefVarError` when an undefined Symbol is found.
"""
function expression_check(allowed_symbols, f_expr; module_name=Utils)
    symbols = parse_expression(f_expr)

    for s in symbols
        if s ∉ allowed_symbols && !isdefined(module_name, s)
            throw(UndefVarError(s))
        end
    end

    return nothing
end

# -----------------------------------------------------------------------------
# Biogeochemistry type construction
# -----------------------------------------------------------------------------

"""
    define_tracer_functions(
        parameters,
        tracers;
        auxiliary_fields=(:PAR,),
        helper_functions=nothing,
        sinking_velocities=nothing,
    )

Create an Oceananigans biogeochemistry model type.

# Arguments
- `parameters`: runtime parameter struct stored as `bgc.parameters`
- `tracers`: `NamedTuple` mapping tracer names (Symbols) to tracer tendency expressions (`Expr`)

# Keywords
- `auxiliary_fields`: tuple of auxiliary field names (Symbols)
- `helper_functions`: optional file path containing helper functions referenced by `tracers`
- `sinking_velocities`: optional `NamedTuple` mapping sinking tracer names to vertical velocity fields
"""
function define_tracer_functions(
    parameters,
    tracers::NamedTuple;
    auxiliary_fields::Tuple=(:PAR,),
    helper_functions=nothing,
    sinking_velocities=nothing,
)
    model_name = gensym(:AgateBGC)

    bgc_type = create_bgc_struct(
        model_name, parameters; sinking_velocities=sinking_velocities
    )

    add_bgc_methods!(
        bgc_type,
        tracers,
        parameters;
        auxiliary_fields=auxiliary_fields,
        helper_functions=helper_functions,
        sinking_velocities=sinking_velocities,
    )

    return bgc_type
end

"""
    create_bgc_struct(struct_name, parameters; sinking_velocities=nothing)

Create a subtype of `AbstractContinuousFormBiogeochemistry` that stores a runtime
parameter struct (and optional sinking velocities).

The generated type is `Adapt.jl`-compatible.

Returns the wrapper type (a `UnionAll`) so it remains valid after `Adapt.adapt`
changes the concrete type parameters.
"""
function create_bgc_struct(struct_name::Symbol, parameters; sinking_velocities=nothing)
    if isnothing(sinking_velocities)
        type_expr = quote
            Base.@kwdef struct $struct_name{PT} <: AbstractContinuousFormBiogeochemistry
                parameters::PT = $parameters
            end
            $struct_name
        end

        T = eval(type_expr)
        eval(:(Adapt.@adapt_structure $struct_name))
        return T
    end

    type_expr = quote
        Base.@kwdef struct $struct_name{PT,W} <: AbstractContinuousFormBiogeochemistry
            parameters::PT = $parameters
            sinking_velocities::W = $sinking_velocities
        end
        $struct_name
    end

    T = eval(type_expr)
    eval(:(Adapt.@adapt_structure $struct_name))
    return T
end

# -----------------------------------------------------------------------------
# Method attachment (public legacy overload + internal implementation)
# -----------------------------------------------------------------------------

@inline function _bgc_wrapper(bgc_type)
    return bgc_type isa UnionAll ? bgc_type : bgc_type.name.wrapper
end

"""
    add_bgc_methods!(
        bgc_type,
        tracers,
        parameters;
        auxiliary_fields=(),
        helper_functions=nothing,
        sinking_velocities=nothing,
    )

Internal implementation: attaches Oceananigans-required methods.

Important:
- Methods are attached to the *wrapper* type so they remain valid after `Adapt.adapt`
  changes type parameters (CPU -> GPU).
- We DO NOT define/override any `OceanBioME.ContinuousBiogeochemistry` trait methods here.
  OceanBioME already provides those; overriding them is version-fragile (e.g., `.model` field
  does not exist in OceanBioME 0.14.0).
"""
function add_bgc_methods!(
    bgc_type,
    tracers::NamedTuple,
    parameters;
    auxiliary_fields::Tuple=(),
    helper_functions=nothing,
    sinking_velocities=nothing,
)
    if !isnothing(helper_functions)
        include(helper_functions)
    end

    coordinates = (:x, :y, :z, :t)
    tracer_vars = keys(tracers)
    aux_field_vars = auxiliary_fields
    all_state_vars = (coordinates..., tracer_vars..., aux_field_vars...)

    wrapper = _bgc_wrapper(bgc_type)

    # Traits on the inner model type (wrapper dispatch survives Adapt.adapt).
    eval(:(required_biogeochemical_tracers(::$(wrapper)) = $(tracer_vars)))
    eval(:(required_biogeochemical_auxiliary_fields(::$(wrapper)) = $(aux_field_vars)))

    # Bind parameter fields into local variables to match tracer expression expectations.
    parameter_fields = fieldnames(typeof(parameters))

    method_vars = Vector{Expr}(undef, length(parameter_fields))
    @inbounds for (i, field) in enumerate(parameter_fields)
        if field in coordinates
            throw(
                ArgumentError("Parameter field name $(field) is reserved for coordinates.")
            )
        end
        method_vars[i] = :($field = bgc.parameters.$field)
    end

    allowed_symbols = (all_state_vars..., parameter_fields...)

    # Define callable tracer tendency methods.
    for (tracer_name, tracer_expression) in pairs(tracers)
        expression_check(allowed_symbols, tracer_expression)

        tracer_method = quote
            function (bgc::$(wrapper))(::Val{$(QuoteNode(tracer_name))}, $(all_state_vars...))
                $(method_vars...)
                return $(tracer_expression)
            end
        end

        eval(tracer_method)
    end

    # Optional sinking velocities.
    if sinking_velocities !== nothing
        sinking_fields = fieldnames(typeof(sinking_velocities))

        # Fallback: no drift unless explicitly provided.
        eval(
            quote
                function biogeochemical_drift_velocity(
                    ::$(wrapper), ::Val{tracer_name}
                ) where {tracer_name}
                    return (u=ZeroField(), v=ZeroField(), w=ZeroField())
                end
            end,
        )

        for tracer_name in sinking_fields
            sinking_method = quote
                function biogeochemical_drift_velocity(
                    bgc::$(wrapper), ::Val{$(QuoteNode(tracer_name))}
                )
                    return (u=ZeroField(), v=ZeroField(), w=bgc.sinking_velocities.$tracer_name)
                end
            end
            eval(sinking_method)
        end
    end

    return bgc_type
end

##############################################################################
# Parameter Utils
###############################################################################

@inline function param_check_length(name::Symbol, expected::Int, got::Int)
    if expected != got
        throw(ArgumentError("$(name) must have length $(expected) but has length $(got)"))
    end
    return nothing
end

@inline function param_check_square_matrix(name::Symbol, n::Int, M::AbstractMatrix)
    if size(M, 1) != n || size(M, 2) != n
        throw(ArgumentError("$(name) must have size ($(n), $(n))"))
    end
    return nothing
end

@inline function param_cast_matrix(::Type{FT}, M::AbstractMatrix) where {FT<:AbstractFloat}
    return eltype(M) === FT ? M : FT.(M)
end

function param_compute_diameters(
    ::Type{FT}, n::Int, spec::DiameterRangeSpecification
) where {FT<:AbstractFloat}
    min_d = FT(spec.min_diameter)
    max_d = FT(spec.max_diameter)

    if n == 1
        return FT[min_d]
    end

    diameters = Vector{FT}(undef, n)

    if spec.splitting === :log_splitting
        log_min = log(min_d)
        log_max = log(max_d)
        step = (log_max - log_min) / FT(n - 1)
        @inbounds for i in 1:n
            diameters[i] = exp(log_min + FT(i - 1) * step)
        end
    elseif spec.splitting === :linear_splitting
        step = (max_d - min_d) / FT(n - 1)
        @inbounds for i in 1:n
            diameters[i] = min_d + FT(i - 1) * step
        end
    else
        throw(ArgumentError("Unsupported splitting method: $(spec.splitting)"))
    end

    return diameters
end

function param_compute_diameters(
    ::Type{FT}, n::Int, spec::DiameterListSpecification
) where {FT<:AbstractFloat}
    _check_length(:diameters, n, length(spec.diameters))
    diameters = Vector{FT}(undef, n)
    @inbounds for i in 1:n
        diameters[i] = FT(spec.diameters[i])
    end
    return diameters
end

end # module
