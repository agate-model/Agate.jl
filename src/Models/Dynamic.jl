"""
Module to dynamically create AbstractContinuousFormBiogeochemistry types (imported from
Oceananigans.Biogeochemistry).
"""

module Dynamic

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

export create_struct, construct_bgc_model


"""
    create_struct(struct_name, priors::Dict) -> AbstractContinuousFormBiogeochemistry

Create a struct of type AbstractContinuousFormBiogeochemistry. Uses field names and default
values defined in priors (which can be, for example, a Dict or NamedTuple).

# Arguments
- `struct_name`: name for the new struct
- `priors`: named sequence of values of the form (field name = default value, ...)
"""
function create_struct(struct_name, priors)

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
    eval(exp)
end


"""
    construct_bgc_model(bgc_struct, tracers, auxiliary_fields=[])

Add required methods to bgc_model struct:
    - method per tracer
    - required_biogeochemical_tracers
    - required_biogeochemical_auxiliary_fields
Optionally can include auxiliary fields.

# Arguments
- `bgc_struct``: struct (of type AbstractContinuousFormBiogeochemistry)
- `tracers`: Dict of named function expressions of the form (name = expression, ...)
- `auxiliary_fields`: optional iterable of auxiliary field variables

# Example
```julia
struct LV <: AbstractContinuousFormBiogeochemistry
    α
    β
    δ
    γ
end
bgc_struct = LV(2/3, 4/3, 1, 1)
tracers = Dict(
    "R" => :(α*R - β*R*F),
    "F" => :(-γ*F + δ*R*F)
)
model = construct_bgc_model(bgc_struct, tracers)
```
"""
function construct_bgc_model(bgc_struct, tracers; auxiliary_fields=[])

    # all BGC models require these base variables
    base_vars = [:x, :y, :z, :t]
    # use collect here in case tracers are NamedTuple rather than Dict
    tracer_vars = Symbol.(collect(keys(tracers)))
    aux_field_vars = Symbol.(auxiliary_fields)
    all_state_vars = vcat(base_vars , tracer_vars, aux_field_vars)

    eval(:(required_biogeochemical_tracers(::$(typeof(bgc_struct))) = Tuple(tracer_vars)))
    eval(:(required_biogeochemical_auxiliary_fields(::$(typeof(bgc_struct))) = Tuple(aux_field_vars)))

    # the variable expressions are evaluated inside the tracer methods built below
    method_vars = []
    for field in fieldnames(typeof(bgc_struct))
        exp = :($field = bgc.$field)
        push!(method_vars, exp)
    end

    for (tracer_name, func_expression) in pairs(tracers)
        tracer_method = quote
            function (bgc::$(typeof(bgc_struct)))(::Val{Symbol($tracer_name)}, $(all_state_vars...))
                $(method_vars...)
                return $(func_expression)
            end
        end
        eval(tracer_method)
    end

    return bgc_struct
end

end #module
