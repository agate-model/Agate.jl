"""
Module to dynamically create Oceananigans.Biogeochmistry types.
"""

module Dynamic

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

include("Library.jl")

export DynamicBGC, construct_bgc_model


struct DynamicBGC <: AbstractContinuousFormBiogeochemistry end


"""
Construct a concrete Oceananigans biogeochemical model instance from a dictionary of
tracer equations. Optionally can also include additonal auxiliary fields.
"""
function construct_bgc_model(tracers; auxiliary_fields=[], priors=Dict())

    for (k,v) in priors
        exp = :($k = $v)
        eval(exp)
    end

    # all BGC models require these base variables
    base_vars = [:x, :y, :z, :t]
    tracer_vars = Symbol.(keys(tracers))
    aux_field_vars = Symbol.(auxiliary_fields)
    all_state_vars = vcat(base_vars , tracer_vars, aux_field_vars)

    required_biogeochemical_tracers(::DynamicBGC) = Tuple(tracer_vars)
    required_biogeochemical_auxiliary_fields(::DynamicBGC) = Tuple(aux_field_vars)

    for (tracer_name, func_expression) in tracers
        tracer_method = quote
            function (::DynamicBGC)(::Val{Symbol($tracer_name)}, $(all_state_vars...))
                return $(func_expression)
            end
        end
        eval(tracer_method)
    end

    return DynamicBGC()
end

end #module
