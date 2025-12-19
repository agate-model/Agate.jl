"""Utilities to construct Oceananigans biogeochemistry types.

Agate constructs `AbstractContinuousFormBiogeochemistry` subtypes from:

- a runtime parameter struct stored as a single field (`bgc.parameters`), and
- construction-time tracer tendency expressions (`Expr`) keyed by tracer name.

The generated types are compatible with Oceananigans and OceanBioME. Runtime parameter
structs can be adapted between CPU and GPU with `Adapt.jl`.
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

export add_bgc_methods!, create_bgc_struct, define_tracer_functions, expression_check

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
- `parameters`: a runtime parameter struct (stored as `bgc.parameters`)
- `tracers`: a `NamedTuple` mapping tracer names (Symbols) to tracer tendency expressions (`Expr`)

# Keywords
- `auxiliary_fields`: a tuple of auxiliary field names (Symbols)
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
    bgc_type = create_bgc_struct(model_name, parameters; sinking_velocities=sinking_velocities)

    add_bgc_methods!(
        bgc_type,
        tracers;
        auxiliary_fields=auxiliary_fields,
        helper_functions=helper_functions,
        include_sinking=sinking_velocities !== nothing,
    )

    return bgc_type
end

"""
    create_bgc_struct(struct_name, parameters; sinking_velocities=nothing)

Create a subtype of `AbstractContinuousFormBiogeochemistry` that stores a runtime
parameter struct (and optional sinking velocities).

The generated type is `Adapt.jl`-compatible.
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

        # --- FIX: return a concrete instantiated type, not the UnionAll. ---
        return T{typeof(parameters)}
    end

    type_expr = quote
        Base.@kwdef struct $struct_name{PT, W} <: AbstractContinuousFormBiogeochemistry
            parameters::PT = $parameters
            sinking_velocities::W = $sinking_velocities
        end
        $struct_name
    end

    T = eval(type_expr)
    eval(:(Adapt.@adapt_structure $struct_name))

    # --- FIX: return a concrete instantiated type, not the UnionAll. ---
    return T{typeof(parameters), typeof(sinking_velocities)}
end

# Compute the index of a field by name without relying on `Base.fieldindex` being in scope.
# This is CPU-side construction logic (not kernel-reachable).
@inline function _field_index(T::Type, name::Symbol)
    names = fieldnames(T)
    @inbounds for i in 1:length(names)
        if names[i] === name
            return i
        end
    end
    throw(ArgumentError("Type $(T) has no field named $(name)."))
end

"""
    add_bgc_methods!(
        bgc_type,
        tracers;
        auxiliary_fields=(),
        helper_functions=nothing,
        include_sinking=false,
    )

Attach Oceananigans-required methods to `bgc_type`.

This function defines:
- `required_biogeochemical_tracers`
- `required_biogeochemical_auxiliary_fields`
- one callable method per tracer tendency
- optionally, `biogeochemical_drift_velocity` for sinking tracers

Each tracer tendency method binds fields of `bgc.parameters` to local variables before
evaluating the provided tracer expression.
"""
function add_bgc_methods!(
    bgc_type,
    tracers::NamedTuple;
    auxiliary_fields::Tuple=(),
    helper_functions=nothing,
    include_sinking::Bool=false,
)
    if !isnothing(helper_functions)
        include(helper_functions)
    end

    coordinates = (:x, :y, :z, :t)
    tracer_vars = keys(tracers)
    aux_field_vars = auxiliary_fields
    all_state_vars = (coordinates..., tracer_vars..., aux_field_vars...)

    eval(:(required_biogeochemical_tracers(::$(bgc_type)) = $(tracer_vars)))
    eval(:(required_biogeochemical_auxiliary_fields(::$(bgc_type)) = $(aux_field_vars)))

    # World-age safe and now concrete: bgc_type has a concrete parameters field type.
    parameters_index = _field_index(bgc_type, :parameters)
    parameter_type = fieldtype(bgc_type, parameters_index)
    parameter_fields = fieldnames(parameter_type)

    method_vars = Vector{Expr}(undef, length(parameter_fields))
    @inbounds for (i, field) in enumerate(parameter_fields)
        if field in coordinates
            throw(ArgumentError("Parameter field name $(field) is reserved for coordinates."))
        end
        method_vars[i] = :($field = bgc.parameters.$field)
    end

    allowed_symbols = (all_state_vars..., parameter_fields...)

    for (tracer_name, tracer_expression) in pairs(tracers)
        expression_check(allowed_symbols, tracer_expression)

        tracer_method = quote
            function (bgc::$(bgc_type))(::Val{$(QuoteNode(tracer_name))}, $(all_state_vars...))
                $(method_vars...)
                return $(tracer_expression)
            end
        end

        eval(tracer_method)
    end

    if include_sinking
        fallback_method = quote
            function biogeochemical_drift_velocity(
                ::$(bgc_type),
                ::Val{tracer_name},
            ) where {tracer_name}
                return (u=ZeroField(), v=ZeroField(), w=ZeroField())
            end
        end
        eval(fallback_method)

        sinking_index = _field_index(bgc_type, :sinking_velocities)
        sinking_type = fieldtype(bgc_type, sinking_index)
        sinking_fields = fieldnames(sinking_type)

        for tracer_name in sinking_fields
            sinking_method = quote
                function biogeochemical_drift_velocity(
                    bgc::$(bgc_type),
                    ::Val{$(QuoteNode(tracer_name))},
                )
                    return (u=ZeroField(), v=ZeroField(), w=bgc.sinking_velocities.$tracer_name)
                end
            end
            eval(sinking_method)
        end
    end

    return bgc_type
end

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

end # module
