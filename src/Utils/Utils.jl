"""Utilities to construct Oceananigans biogeochemistry types.

Agate builds `AbstractContinuousFormBiogeochemistry` subtypes from:

- a runtime `parameters` struct (stored as a single field), and
- tracer tendency expressions (construction-time `Expr` values).

The generated biogeochemistry types are compatible with Oceananigans and OceanBioME.
"""

module Utils

using Agate.Library.Mortality
using Agate.Library.Nutrients
using Agate.Library.Photosynthesis
using Agate.Library.Predation
using Agate.Library.Remineralization

using UUIDs

using OceanBioME
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Fields: ZeroField

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

export define_tracer_functions, expression_check, create_bgc_struct, add_bgc_methods!

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
- `parameters`: a runtime parameter struct (stored as a single field of the generated type)
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
    model_name = Symbol(uuid1())
    bgc_model = create_bgc_struct(model_name, parameters; sinking_velocities=sinking_velocities)

    add_bgc_methods!(
        bgc_model,
        tracers;
        auxiliary_fields=auxiliary_fields,
        helper_functions=helper_functions,
        include_sinking=sinking_velocities !== nothing,
    )

    return bgc_model
end

"""
    create_bgc_struct(struct_name, parameters; sinking_velocities=nothing)

Create a subtype of `AbstractContinuousFormBiogeochemistry` that stores a runtime
parameter struct (and optional sinking velocities).

The generated type is returned.
"""
function create_bgc_struct(struct_name::Symbol, parameters; sinking_velocities=nothing)
    if isnothing(sinking_velocities)
        exp = quote
            Base.@kwdef struct $struct_name{PT} <: AbstractContinuousFormBiogeochemistry
                parameters::PT = $parameters
            end
            $struct_name
        end
        return eval(exp)
    else
        exp = quote
            Base.@kwdef struct $struct_name{PT, W} <: AbstractContinuousFormBiogeochemistry
                parameters::PT = $parameters
                sinking_velocities::W = $sinking_velocities
            end
            $struct_name
        end
        return eval(exp)
    end
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

The generated tracer methods bind each field of `bgc.parameters` to a local variable
before evaluating the provided tracer expression.
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
    tracer_vars = Tuple(Symbol.(keys(tracers)))
    aux_field_vars = Tuple(Symbol.(auxiliary_fields))

    all_state_vars = (coordinates..., tracer_vars..., aux_field_vars...)

    eval(:(required_biogeochemical_tracers(::$(bgc_type)) = $(tracer_vars,)))
    eval(:(required_biogeochemical_auxiliary_fields(::$(bgc_type)) = $(aux_field_vars,)))

    defaults = bgc_type()
    parameter_type = typeof(defaults.parameters)
    parameter_fields = fieldnames(parameter_type)

    method_vars = Expr[]
    for field in parameter_fields
        if field in coordinates
            throw(ArgumentError("Parameter field name $(field) is reserved for coordinates."))
        end
        push!(method_vars, :($field = bgc.parameters.$field))
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
        sinking_keys = Tuple(keys(defaults.sinking_velocities))

        for tracer in tracer_vars
            if tracer in sinking_keys
                sink_velocity_method = quote
                    function biogeochemical_drift_velocity(
                        bgc::$(bgc_type),
                        ::Val{$(QuoteNode(tracer))},
                    )
                        return (u=ZeroField(), v=ZeroField(), w=bgc.sinking_velocities.$tracer)
                    end
                end
            else
                sink_velocity_method = quote
                    function biogeochemical_drift_velocity(
                        bgc::$(bgc_type),
                        ::Val{$(QuoteNode(tracer))},
                    )
                        return (u=ZeroField(), v=ZeroField(), w=ZeroField())
                    end
                end
            end

            eval(sink_velocity_method)
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
