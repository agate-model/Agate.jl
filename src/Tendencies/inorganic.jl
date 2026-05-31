"""Inorganic pool tendency using the configured growth and nutrient couplings."""
function inorganic_tendency(
    config::TendencyConfig{Formulation,Zooplankton,Limitation};
    target::Symbol,
    remineralization=nothing,
    stoichiometry=nothing,
) where {Formulation,Zooplankton,Limitation}
    nutrients = config.nutrients
    target_nutrient = target in map(tracer_name, nutrients) ? target_coupling(nutrients, target) : nothing
    sources = remineralization === nothing ? remineralization_sources(target_nutrient) : remineralization
    ratio_name = stoichiometry === nothing && target_nutrient !== nothing ? stoichiometry_name(target_nutrient) : something(stoichiometry, :one)

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        resources = nutrient_resources(tracer_values, nutrients)
        half_saturations = half_saturation_parameters(parameters, nutrients)
        uptake = growth_sum(
            Val(Formulation),
            Val(Limitation),
            tracer_values.plankton,
            resources,
            tracer_values.PAR,
            parameters.maximum_growth_rate,
            half_saturations,
            growth_parameters(Val(Formulation), parameters)...,
        )

        remin = remineralization_sum(tracer_values, parameters, sources)
        ratio = stoichiometry_multiplier(parameters, Val(ratio_name))

        return remin - ratio * uptake
    end

    return CompiledEquation(f)
end
