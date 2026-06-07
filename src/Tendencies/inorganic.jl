function inorganic_target_coupling(config, target::Symbol, remineralization, stoichiometry)
    nutrients = config.nutrients
    target_is_nutrient = target in map(tracer_name, nutrients)
    target_nutrient = target_is_nutrient ? target_coupling(nutrients, target) : nothing

    if !target_is_nutrient
        remineralization === nothing &&
            error(
                "Inorganic target $target is not configured as a nutrient; " *
                "pass explicit remineralization sources."
            )
        @warn "Inorganic target is not configured as a nutrient; using explicit inorganic-pool coupling." target=target remineralization=remineralization stoichiometry=stoichiometry
    end

    sources = remineralization === nothing ? remineralization_sources(target_nutrient) : remineralization
    ratio_name = stoichiometry === nothing && target_nutrient !== nothing ?
        stoichiometry_name(target_nutrient) :
        something(stoichiometry, :one)

    return sources, ratio_name
end

"""Inorganic pool tendency with simple-detritus nutrient cycling."""
function inorganic_tendency(
    config::TendencyConfig{Growth,:simple_detritus,Zooplankton,Limitation};
    target::Symbol,
    remineralization=nothing,
    stoichiometry=nothing,
    export_fraction::Symbol=:mortality_export_fraction,
) where {Growth,Zooplankton,Limitation}
    nutrients = config.nutrients
    sources, ratio_name = inorganic_target_coupling(
        config, target, remineralization, stoichiometry
    )

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        resources = nutrient_resources(tracer_values, nutrients)
        half_saturations = half_saturation_parameters(parameters, nutrients)
        uptake = growth_sum(
            Val(Growth),
            Val(Limitation),
            tracer_values.plankton,
            resources,
            tracer_values.PAR,
            parameters.maximum_growth_rate,
            half_saturations,
            growth_parameters(Val(Growth), parameters)...,
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

"""Inorganic pool tendency with DOM/POM nutrient cycling."""
function inorganic_tendency(
    config::TendencyConfig{Growth,:dom_pom,Zooplankton,Limitation};
    target::Symbol,
    remineralization=nothing,
    stoichiometry=nothing,
) where {Growth,Zooplankton,Limitation}
    nutrients = config.nutrients
    sources, ratio_name = inorganic_target_coupling(
        config, target, remineralization, stoichiometry
    )

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        resources = nutrient_resources(tracer_values, nutrients)
        half_saturations = half_saturation_parameters(parameters, nutrients)
        uptake = growth_sum(
            Val(Growth),
            Val(Limitation),
            tracer_values.plankton,
            resources,
            tracer_values.PAR,
            parameters.maximum_growth_rate,
            half_saturations,
            growth_parameters(Val(Growth), parameters)...,
        )

        remin = remineralization_sum(tracer_values, parameters, sources)
        ratio = stoichiometry_multiplier(parameters, Val(ratio_name))

        return remin - ratio * uptake
    end

    return CompiledEquation(f)
end
