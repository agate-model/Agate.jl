"""
Module to dynamically created Oceananigans.Biogeochmistry types from model equations....
"""

module Dynamic

export DynamicBGC, create_bgc_model

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry


"""

"""
struct DynamicBGC <: AbstractContinuousFormBiogeochemistry end


"""

"""
# TODO: add types
function create_bgc_model(tracers, aux_field_vars)
    
    # all BGC models require these base variables
    base_vars = [:x, :y, :z, :t]
    tracer_vars = Symbol.(keys(tracers))
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