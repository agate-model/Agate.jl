"""
Module to dynamically create Oceananigans.Biogeochemistry types.
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
        auxiliary_fields=[:PAR],
        helper_functions=nothing,
        sinking_velocities=nothing,
    )

Create an Oceananigans.Biogeochemistry model type.

# Arguments
- `parameters`: NamedTuple of values of the form ((<field name> = <default value>, ...)
- `tracers`: Dictionary of the form (<tracer name> => <tracer method expression>, ...)

# Keywords
- `auxiliary_fields`: an iterable of auxiliary field variables, defaults to `[:PAR]`
- `helper_functions`: optional path to a file of helper functions used in tracer expressions
- `sinking_velocities`: optional NamedTuple of constant sinking, of fields (i.e. `ZFaceField(...)`)
   for any tracers which sink returned by OceanBioME.Sediments: sinking_tracers.

!!! warning

    Note that the field names defined in `parameters` can't be any of [:x, :y, :z, :t], as these
    are reserved for coordinates, and they must include all parameters used in the `tracers`
    expressions.

!!! warning

    The tracer expressions must use methods that are either defined within this module or
    passed in the `helper_functions` file.

# Example
```julia
using Agate

parameters = (α=2 / 3, β=4 / 3, δ=1, γ=1)
tracers = Dict("R" => :(α * R - β * R * F), "F" => :(-γ * F + δ * R * F))
LV = define_tracer_functions(parameters, tracers)
```
"""
function define_tracer_functions(
    parameters,
    tracers;
    auxiliary_fields=[:PAR],
    helper_functions=nothing,
    sinking_velocities=nothing,
)
    # create a universaly unique identifier (UUID) for the model struct
    model_name = Symbol(uuid1())
    bgc_model = create_bgc_struct(model_name, parameters, sinking_velocities)
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
    create_bgc_struct(
        struct_name,
        parameters;
        sinking_velocities=nothing,
    )

Create a subtype of Oceananigans.Biogeochemistry with field names defined in `parameters`.

# Arguments
- `struct_name`: name for the struct to create passed as a Symbol (the new struct will be
   accessible as `Agate.Models.Biogeochemistry.<struct_name>`)
- `parameters`: NamedTuple of values of the form (<field name> = <default value>, ...)

# Keywords
- `sinking_velocities`: optional NamedTuple of constant sinking, of fields (i.e. `ZFaceField(...)`)
   for any tracers which sink returned by OceanBioME.Sediments: sinking_tracers.

!!! warning

    Note that the field names defined in `parameters` can't be any of [:x, :y, :z, :t] as these
    are reserved for coordinates.

# Example
```julia
using Agate.Models.Biogeochemistry: create_bgc_struct

create_bgc_struct(:LV, (α=2/3, β=4/3,  δ=1, γ=1))
```
"""
function create_bgc_struct(struct_name, parameters, sinking_velocities=nothing)
    # have to create an expression for each struct field
    # this is of the form `<field name>:: <field type> = <field value>`
    # create and store expressions in an array before struct contsruction
    fields = []
    # need to also keep track of all parameter types to return a parametric struct
    type_names = Set()
    for (k, v) in pairs(parameters)
        if k in [:x, :y, :z, :t]
            throw(
                DomainError(k, "field names in parameters can't be any of [:x, :y, :z, :t]")
            )
        end

        type = typeof(v)
        # should we handle ints/floats separately?
        if type <: Real
            type_symbol = :FT
        elseif type <: AbstractVector
            type_symbol = :VT
        elseif type <: AbstractMatrix
            type_symbol = :MT
        else
            error("Unsupported type for field $k")
        end
        push!(type_names, type_symbol)

        exp = Expr(:(=), Expr(:(::), k, type_symbol), v)
        push!(fields, exp)
    end

    # optionally add a field for sinking velocities
    if !isnothing(sinking_velocities)
        # using W here for consistency with OceanBioME
        type_symbol = :W
        push!(type_names, type_symbol)
        exp = Expr(:(=), Expr(:(::), :sinking_velocities, type_symbol), sinking_velocities)
        push!(fields, exp)
    end

    # construct struct (include default parameter values so have to use kwdef)
    exp = quote
        Base.@kwdef struct $struct_name{$(type_names...)} <:
                           AbstractContinuousFormBiogeochemistry
            $(fields...)
        end
        $struct_name # return the type
    end
    return eval(exp)
end

"""
    add_bgc_methods!(
        bgc_type,
        tracers;
        auxiliary_fields=[],
        helper_functions=nothing,
        include_sinking=false,
    )

Add methods to `bgc_type` required of Oceananigans.Biogeochemistry:
    - `required_biogeochemical_tracers`
    - `required_biogeochemical_auxiliary_fields`
    - a method per tracer specifying how it evolves in time
    - optionally adds `biogeochemical_drift_velocity` (if `include_sinking` is true)

!!! info

    before passing the `bgc_type` to Oceananigans, it needs to be wrapped in an
    `OceanBioME.Biogeochemistry()` object, which adds additional methods not defined here:
    - `biogeochemical_auxiliary_fields`
    - `update_biogeochemical_state!`

# Arguments
- `bgc_type`: subtype of Oceananigans.Biogeochemistry (returned by `create_bgc_struct`)
- `tracers`: Dictionary of the form (<tracer name> => <tracer method expression>, ...)

# Keywords
- `auxiliary_fields`: an optional iterable of auxiliary field variables, defaults to `[]`
- `helper_functions`: optional path to a file of helper functions used in tracer expressions
- `include_sinking`: boolean indicating whether the model includes sinking tracers, if true
   adds corresponding OceanBioME methods (e.g., `biogeochemical_drift_velocity()`), defaults
   to false

!!! warning

    Note that the field names of `bgc_type` can't be any of [:x, :y, :z, :t] (as these are reserved
    for coordinates) and they must include all parameters used in the `tracers` expressions.

!!! warning

    The tracer expressions must use methods that are either defined within this module or
    passed in the `helper_functions` file.

# Example
```julia
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using using Agate.Models.Biogeochemistry: add_bgc_methods!

struct LV <: AbstractContinuousFormBiogeochemistry
    α
    β
    δ
    γ
end

tracers = Dict(
    "R" => :(α*R - β*R*F),
    "F" => :(-γ*F + δ*R*F)
)

add_bgc_methods!(LV, tracers)
```
"""
function add_bgc_methods!(
    bgc_type, tracers; auxiliary_fields=[], helper_functions=nothing, include_sinking=false
)
    if !isnothing(helper_functions)
        include(helper_functions)
    end

    coordinates = [:x, :y, :z, :t]
    # use collect here in case tracers are NamedTuple rather than Dict
    tracer_vars = Symbol.(collect(keys(tracers)))
    aux_field_vars = Symbol.(auxiliary_fields)
    all_state_vars = vcat(coordinates, tracer_vars, aux_field_vars)

    eval(:(required_biogeochemical_tracers(::$(bgc_type)) = $(tracer_vars...,)))
    eval(:(required_biogeochemical_auxiliary_fields(::$(bgc_type)) = $(aux_field_vars...,)))

    # the BGC struct holds all the user defined parameters in its fields
    params = fieldnames(bgc_type)
    # the tracer methods rely on model parameters to be declared as variables in the method
    # scope, here we create expressions that define these local variables
    method_vars = []
    for param in params
        if param in [:x, :y, :z, :t]
            throw(
                DomainError(field, "$bgc_type field names can't be any of [:x, :y, :z, :t]")
            )
        end
        # the expressions are evaluated inside the tracer methods below, which take in a
        # `bgc` object (of bgc_type)
        exp = :($param = bgc.$param)
        push!(method_vars, exp)
    end

    # create core tracer methods
    for (tracer_name, tracer_expression) in pairs(tracers)

        # throws an exception if there are any issues with tracer_expression (see docstring)
        expression_check(vcat(all_state_vars, collect(params)), tracer_expression)

        tracer_method = quote
            function (bgc::$(bgc_type))(::Val{Symbol($tracer_name)}, $(all_state_vars...))
                $(method_vars...)
                return $(tracer_expression)
            end
        end
        eval(tracer_method)
    end

    # set up tracer sinking methods
    if include_sinking
        # `biogeochemical_drift_velocity` is an optional Oceananigans.Biogeochemistry method
        # that returns a NamedTuple of velocity fields for a tracer with keys `u`, `v`, `w`
        # implementation here based on:
        # https://github.com/OceanBioME/OceanBioME.jl/blob/a84ab98465b220e805232c46015759a6d5280536/src/Models/AdvectedPopulations/NPZD.jl#L286
        sink_velocity_method = quote
            function biogeochemical_drift_velocity(
                bgc::$(bgc_type), ::Val{tracer_name}
            ) where {tracer_name}
                if tracer_name in keys(bgc.sinking_velocities)
                    return (u=ZeroField(), v=ZeroField(), w=bgc.sinking_velocities[tracer_name])
                else
                    return (u=ZeroField(), v=ZeroField(), w=ZeroField())
                end
            end
        end
        eval(sink_velocity_method)
    end

    return bgc_type
end

"""
    parse_expression(f_expr) -> Vector

Return all symbols (function names and argument names) called in expression.

!!! info

    The input here is expected to be a tracer method defined inside a quote (see example).
    This allows us to check whether all functions and variables in the tracer method are defined.

# Example
```julia
using Agate.Utils: parse_expression

parse_expression(:(α * x - β * x * y))
```
"""
function parse_expression(f_expr)
    symbols = []

    expressions = [f_expr]
    for exp in expressions
        # if the arg isn't a symbol or another expression, it's a value --> ignore
        for arg in exp.args
            if isa(arg, Expr)
                push!(expressions, arg)
            elseif isa(arg, Symbol)
                push!(symbols, arg)
            end
        end
    end

    return symbols
end

"""
    expression_check(args, f_expr) -> nothing

Check that all methods and arguments are defined. Specifically:
    - vector `args` contains all arguments of expression `f_expr`
    - all methods called in `f_expr` are defined in module (e.g., Base, Main, Agate)
If not, throws an UnderVarError.
"""
function expression_check(args, f_expr; module_name=Utils)
    symbols = parse_expression(f_expr)
    for s in symbols
        if s ∉ args && !isdefined(module_name, s)
            throw(UndefVarError(s))
        end
    end
end

end #module
