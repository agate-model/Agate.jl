"""
Module to dynamically create Oceananigans.Biogeochmistry types.
"""

module Dynamic

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

# all Agate functions have to be in scope to evaluate expressions that might contain them
include("Library.jl")

export create_struct, construct_bgc_model


"""
Create struct of type AbstractContinuousFormBiogeochemistry...
"""
function create_struct(priors, struct_name)

    fields = []
    for (k,v) in priors
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
Construct a concrete Oceananigans biogeochemical model instance from a dictionary of
tracer equations. Optionally can include auxiliary fields...
"""
function construct_bgc_model(tracers, bgc_struct; auxiliary_fields=[])

    # all BGC models require these base variables
    base_vars = [:x, :y, :z, :t]
    tracer_vars = Symbol.(keys(tracers))
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

    for (tracer_name, func_expression) in tracers
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
