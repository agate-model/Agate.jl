"""Inorganic pool tendency for Smith-detritus nutrient cycling."""
function inorganic_tendency(
    config::TendencyConfig{:smith_detritus,Zooplankton,Limitation};
    target::Symbol,
    remineralization=nothing,
    stoichiometry=nothing,
    export_fraction::Symbol=:mortality_export_fraction,
) where {Zooplankton,Limitation}
    nutrients = config.nutrients
    target_nutrient = target in map(tracer_name, nutrients) ? target_coupling(nutrients, target) : nothing
    sources = remineralization === nothing ? remineralization_sources(target_nutrient) : remineralization
    ratio_name = stoichiometry === nothing && target_nutrient !== nothing ? stoichiometry_name(target_nutrient) : something(stoichiometry, :one)

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        resources = nutrient_resources(tracer_values, nutrients)
        half_saturations = half_saturation_parameters(parameters, nutrients)
        uptake = growth_sum(
            Val(:smith_detritus),
            Val(Limitation),
            tracer_values.plankton,
            resources,
            tracer_values.PAR,
            parameters.maximum_growth_rate,
            half_saturations,
            growth_parameters(Val(:smith_detritus), parameters)...,
        )

        mortality = mortality_loss_sum(
            tracer_values.plankton,
            parameters.linear_mortality,
            parameters.quadratic_mortality,
        )
        export_frac = parameter_value(parameters, Val(export_fraction))
        exported_mortality = export_frac * mortality

        remin = remineralization_sum(tracer_values, parameters, sources)
        ratio = stoichiometry_multiplier(parameters, Val(ratio_name))

        return exported_mortality + remin - ratio * uptake
    end

    return CompiledEquation(f)
end

"""Inorganic pool tendency for Geider DOM/POM nutrient cycling."""
function inorganic_tendency(
    config::TendencyConfig{:geider_dom_pom,Zooplankton,Limitation};
    target::Symbol,
    remineralization=nothing,
    stoichiometry=nothing,
) where {Zooplankton,Limitation}
    nutrients = config.nutrients
    target_nutrient = target in map(tracer_name, nutrients) ? target_coupling(nutrients, target) : nothing
    sources = remineralization === nothing ? remineralization_sources(target_nutrient) : remineralization
    ratio_name = stoichiometry === nothing && target_nutrient !== nothing ? stoichiometry_name(target_nutrient) : something(stoichiometry, :one)

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        resources = nutrient_resources(tracer_values, nutrients)
        half_saturations = half_saturation_parameters(parameters, nutrients)
        uptake = growth_sum(
            Val(:geider_dom_pom),
            Val(Limitation),
            tracer_values.plankton,
            resources,
            tracer_values.PAR,
            parameters.maximum_growth_rate,
            half_saturations,
            growth_parameters(Val(:geider_dom_pom), parameters)...,
        )

        remin = remineralization_sum(tracer_values, parameters, sources)
        ratio = stoichiometry_multiplier(parameters, Val(ratio_name))

        return remin - ratio * uptake
    end

    return CompiledEquation(f)
end
