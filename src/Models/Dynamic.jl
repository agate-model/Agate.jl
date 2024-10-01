"""
Module to dynamically create AbstractContinuousFormBiogeochemistry types (imported from
Oceananigans.Biogeochemistry).
"""

module Dynamic

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields

export create_bgc_struct, add_bgc_methods


"""
    create_bgc_type(struct_name, priors::Dict) -> DataType

Create a subtype of AbstractContinuousFormBiogeochemistry. Uses field names and default
values defined in priors (which can be, for example, a Dict or NamedTuple).

# Arguments
- `struct_name`: name for the new struct
- `priors`: named sequence of values of the form (field name = default value, ...)

Note that the field names in priors can't be any of [:x,:y,:z,:t].

# Example
create_bgc_struct(:LV, (α=2/3, β=4/3,  δ=1, γ=1))
"""
function create_bgc_struct(struct_name, priors)

    fields = []
    for (k,v) in pairs(priors)
        if k in [:x,:y,:z,:t]
            throw(DomainError(k, "field names in priors can't be any of [:x,:y,:z,:t]"))
        end
        exp = :($k = $v)
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
    add_bgc_methods(bgc_type, tracers, auxiliary_fields=[], helper_functions=()) -> DataType

Add methods to bgc_type required of AbstractContinuousFormBiogeochemistry:
    - required_biogeochemical_tracers
    - required_biogeochemical_auxiliary_fields
    - a method per tracer

# Arguments
- `bgc_type`: subtype of AbstractContinuousFormBiogeochemistry
- `tracers`: dictionary of the form (name => expression, ...)
- `auxiliary_fields`: optional iterable of auxiliary field variables
- `helper_functions`: optional path to a file of helper functions used in tracer expressions

Note that the field names of bgc_type can't be any of [:x,:y,:z,:t] and they must include
all parameters used in the tracers expressions. The expressions must use methods that are
either defined within this module or passed in the helper_functions file.

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

add_bgc_methods(LV, tracers)
```
"""
function add_bgc_methods(bgc_type, tracers; auxiliary_fields=[], helper_functions="")

    if helper_functions != ""
        include(helper_functions)
    end

    # all BGC models require these 4 base variables
    base_vars = [:x, :y, :z, :t]
    # use collect here in case tracers are NamedTuple rather than Dict
    tracer_vars = Symbol.(collect(keys(tracers)))
    println(tracer_vars)
    aux_field_vars = Symbol.(auxiliary_fields)
    all_state_vars = vcat(base_vars , tracer_vars, aux_field_vars)
    println()

    eval(:(required_biogeochemical_tracers(::$(bgc_type)) = $tracer_vars))
    eval(:(required_biogeochemical_auxiliary_fields(::$(bgc_type)) = $aux_field_vars))

    params = fieldnames(bgc_type)
    method_vars = []
    for param in params
        if param in [:x,:y,:z,:t]
            throw(DomainError(field, "$bgc_type field names can't be any of [:x,:y,:z,:t]"))
        end
        # the expressions are evaluated inside the tracer methods below, which take in a
        # `bgc` object (of bgc_type)
        exp = :($param = bgc.$param)
        push!(method_vars, exp)
    end

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

    # Q: this return statement is unnecessary - is it good practice to leave or remove it?
    return bgc_type
end


"""
    parse_expression(f_expr) -> Tuple(Vector, Vector)

Returns all argument names and methods called in expression.

# Example
```Julia
parse_expression(:(α * x - β * x * y))
"""
function parse_expression(f_expr)
    argnames = []
    methods = []

    expressions = [f_expr]
    for exp in expressions
        args = exp.args
        push!(methods, args[1])
        # if the arg isn't a symbol or another expression, it's a value --> ignore
        for arg in args[2:end]
            if typeof(arg) == Expr
                push!(expressions, arg)
            elseif typeof(arg) == Symbol
                push!(argnames, arg)
            end
        end
    end

    return argnames, methods
end


"""
    expression_check(params, f_expr) -> nothing

Checks that all methods and arguments are defined. Specifically:
    - vector params contains all arguments of expression f_expr
    - all methods called in expression are defined in module (e.g., Base, Main, Agate)
If not, throws an UnderVarError.
"""
function expression_check(params, f_expr; module_name=Dynamic)
    argnames, methods = parse_expression(f_expr)
    for arg in argnames
        if arg ∉ params
            throw(UndefVarError(arg))
        end
    end
    for m in methods
        if !isdefined(module_name, m)
            throw(UndefVarError(m))
        end
    end
end

end #module
