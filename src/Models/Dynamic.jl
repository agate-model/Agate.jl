"""
Module to dynamically create Oceananigans.Biogeochemistry types.
"""

module Dynamic

using UUIDs

using OceanBioME: setup_velocity_fields

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Fields: ZeroField

import Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity
import OceanBioME.Models.Sediments: sinking_tracers

export define_tracer_functions

"""
    define_tracer_functions(
        parameters,
        tracers;
        auxiliary_fields=[:PAR],
        helper_functions=nothing,
        sinking_tracers=nothing,
        grid=nothing,
        open_bottom=false,
    ) -> DataType

Creates an Oceananigans.Biogeochemistry model type.

# Arguments
- `parameters`: named sequence of values of the form ((<field name> = <default value>, ...)
- `tracers`: dictionary of the form (<name> => <expression>, ...)

# Keywords
- `auxiliary_fields`: an iterable of auxiliary field variables, defaults to `[:PAR,]`
- `helper_functions`: optional path to a file of helper functions used in tracer expressions
- `sinking_tracers`: optional NamedTuple of sinking speeds (passed as positive values) of
   the form (<tracer name> = <speed>, ...)
- `grid`: optional Oceananigans grid object defining the geometry to build the model on, must
   be passed if `sinking_tracers` is defined
- `open_bottom`: indicates whether the sinking velocity should be smoothly brought to zero
   at the bottom to prevent the tracers leaving the domain, defaults to `true`, which means
   the bottom is open and the tracers leave (i.e., no slowing of velocity to 0 is applied)

Note that the field names defined in `parameters` can't be any of [:x, :y, :z, :t] (as these
are reserved for coordinates) and they must include all parameters used in the `tracers`
expressions. The expressions must use methods that are either defined within this module or
passed in the `helper_functions` file.

# Example
```julia
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
    sinking_tracers=nothing,
    grid=nothing,
    open_bottom=true,
)
    # create a universaly unique identifier (UUID) for the model struct
    model_name = Symbol(uuid1())
    bgc_model = create_bgc_struct(
        model_name, parameters, sinking_tracers, grid, open_bottom
    )
    add_bgc_methods!(
        bgc_model,
        tracers;
        auxiliary_fields=auxiliary_fields,
        helper_functions=helper_functions,
        sinking_tracers=sinking_tracers,
    )
    return bgc_model
end

"""
    create_bgc_struct(
        struct_name,
        parameters;
        sinking_tracers=nothing,
        grid=nothing,
        open_bottom=false,
    ) -> DataType

Create a subtype of Oceananigans.Biogeochemistry with field names defined in `parameters`.

# Arguments
- `struct_name`: name for the new struct passed as a Symbol. The struct will be accessible
   as: `Agate.Models.Dynamic.<struct_name>`
- `parameters`: named sequence of values of the form (<field name> = <default value>, ...)

# Keywords
- `sinking_tracers`: optional NamedTuple of sinking speeds (passed as positive values) of
   the form (<tracer name> = <speed>, ...)
- `grid`: optional Oceananigans grid object defining the geometry to build the model on, must
   be passed if `sinking_tracers` is defined
- `open_bottom`: indicates whether the sinking velocity should be smoothly brought to zero
   at the bottom to prevent the tracers leaving the domain, defaults to `true`, which means
   the bottom is open and the tracers leave (i.e., no slowing of velocity to 0 is applied)

Note that the field names defined in `parameters` can't be any of [:x, :y, :z, :t] as these
are reserved for coordinates.

# Example
```julia
create_bgc_struct(:LV, (α=2/3, β=4/3,  δ=1, γ=1))
```
"""
function create_bgc_struct(
    struct_name, parameters, sinking_tracers=nothing, grid=nothing, open_bottom=nothing
)
    fields = []
    for (k, v) in pairs(parameters)
        if k in [:x, :y, :z, :t]
            throw(
                DomainError(k, "field names in parameters can't be any of [:x, :y, :z, :t]")
            )
        end
        exp = :($k = $v)
        push!(fields, exp)
    end

    # additional parameters required if tracers sink
    if !isnothing(sinking_tracers)
        if isnothing(grid)
            throw(ArgumentError("grid must be defined to setup tracer sinking"))
        end
        # `setup_velocity_fields` multiplies tracer speeds by -1 when creating the fields
        sinking_velocities = setup_velocity_fields(sinking_tracers, grid, open_bottom)
        exp = :(sinking_velocities = $(sinking_velocities))
        push!(fields, exp)
    end

    exp = quote
        Base.@kwdef struct $struct_name <: AbstractContinuousFormBiogeochemistry
            $(fields...)
        end
    end
    return eval(exp)
end

"""
    add_bgc_methods!(
        bgc_type,
        tracers;
        auxiliary_fields=[:PAR],
        helper_functions=nothing,
        sinking_tracers=nothing,
    ) -> DataType

Add methods to `bgc_type` required of `AbstractContinuousFormBiogeochemistry`:
    - `required_biogeochemical_tracers`
    - `required_biogeochemical_auxiliary_fields`
    - a method per tracer specifying how it evolves in time
    - optionally adds `biogeochemical_drift_velocity` (if `sinking_tracers` is defined)

WARNING: `biogeochenical_auxiliary_fields` method must also be defined to make use of the
`auxiliary_fields`. This method is added when `OceanBioME.Biogeochemistry(bgc_type())` is
instantiated alongside an `update_biogeochemical_state!` method.

# Arguments
- `bgc_type`: subtype of `AbstractContinuousFormBiogeochemistry` (returned by `create_bgc_struct`)
- `tracers`: dictionary of the form (<name> => <expression>, ...)

# Keywords
- `auxiliary_fields`: an optional iterable of auxiliary field variables
- `helper_functions`: optional path to a file of helper functions used in tracer expressions
- `sinking_tracers`: optional NamedTuple of sinking speeds of the form
   (<tracer name> = <speed>, ...), convention is that speeds are defined as positive values

Note that the field names of `bgc_type` can't be any of [:x, :y, :z, :t] (as these are reserved
for coordinates) and they must include all parameters used in the `tracers` expressions. The
expressions must use methods that are either defined within this module or passed in the
`helper_functions` file.

# Example
```julia
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

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
    bgc_type,
    tracers;
    auxiliary_fields=[],
    helper_functions=nothing,
    sinking_tracers=nothing,
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

    params = fieldnames(bgc_type)
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

    # set up tracer sinking
    if !isnothing(sinking_tracers)

        # `biogeochemical_drift_velocity` is an optional Oceananigans.Biogeochemistry method
        # that returns a NamedTuple of velocity fields for a tracer with keys `u`, `v`, `w`
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

        # this function is used in the OceanBioME sediment models
        eval(:(sinking_tracers(bgc::$(bgc_type)) = keys(bgc.sinking_velocities)))
    end

    return bgc_type
end

"""
    parse_expression(f_expr) -> Vector

Returns all symbols (argument names and method names) called in expression.

# Example
```julia
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

Checks that all methods and arguments are defined. Specifically:
    - vector `args` contains all arguments of expression `f_expr`
    - all methods called in `f_expr` are defined in module (e.g., Base, Main, Agate)
If not, throws an UnderVarError.
"""
function expression_check(args, f_expr; module_name=Dynamic)
    symbols = parse_expression(f_expr)
    for s in symbols
        if s ∉ args && !isdefined(module_name, s)
            throw(UndefVarError(s))
        end
    end
end

end #module
