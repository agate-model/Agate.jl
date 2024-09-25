"""
Module to dynamically create AbstractContinuousFormBiogeochemistry types (imported from
Oceananigans.Biogeochemistry).
"""

module Dynamic

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

export create_bgc_struct, add_bgc_methods


"""
    create_bgc_type(struct_name, priors::Dict) -> DataType

Create a subtype of AbstractContinuousFormBiogeochemistry. Uses field names and default
values defined in priors (which can be, for example, a Dict or NamedTuple).

# Arguments
- `struct_name`: name for the new struct
- `priors`: named sequence of values of the form (field name = default value, ...)

# Example
create_bgc_struct(:LV, (α=2/3, β=4/3,  δ=1, γ=1))
"""
function create_bgc_struct(struct_name, priors)

    fields = []
    for (k,v) in pairs(priors)
        exp = :($k = $v)
        push!(fields, exp)
    end

    exp = quote
        @kwdef struct $struct_name <: AbstractContinuousFormBiogeochemistry
            $(fields...)
        end
    end
    return eval(exp)
end


"""
    add_bgc_methods(bgc_type, tracers, auxiliary_fields=[]) -> DataType

Add required methods to bgc_type:
    - method per tracer
    - required_biogeochemical_tracers
    - required_biogeochemical_auxiliary_fields

# Arguments
- `bgc_type`: subtype of AbstractContinuousFormBiogeochemistry
- `tracers`: Dict of named function expressions of the form (name = expression, ...)
- `auxiliary_fields`: optional iterable of auxiliary field variables

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
function add_bgc_methods(bgc_type, tracers; auxiliary_fields=[])

    # all BGC models require these base variables
    base_vars = [:x, :y, :z, :t]
    # use collect here in case tracers are NamedTuple rather than Dict
    tracer_vars = Symbol.(collect(keys(tracers)))
    aux_field_vars = Symbol.(auxiliary_fields)
    all_state_vars = vcat(base_vars , tracer_vars, aux_field_vars)

    eval(:(required_biogeochemical_tracers(::$(bgc_type)) = Tuple(tracer_vars)))
    eval(:(required_biogeochemical_auxiliary_fields(::$(bgc_type)) = Tuple(aux_field_vars)))

    # the variable expressions are evaluated inside the tracer methods built below, which
    # take in a `bgc` object (of bgc_type)
    method_vars = []
    for field in fieldnames(bgc_type)
        exp = :($field = bgc.$field)
        push!(method_vars, exp)
    end

    for (tracer_name, func_expression) in pairs(tracers)
        tracer_method = quote
            function (bgc::$(bgc_type))(::Val{Symbol($tracer_name)}, $(all_state_vars...))
                $(method_vars...)
                return $(func_expression)
            end
        end
        eval(tracer_method)
    end

    # Q: this return statement is unnecessary - is it good practice to leave it or remove it?
    return bgc_type
end

end #module
