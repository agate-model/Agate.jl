"""Phytoplankton tendency configured by growth and nutrient limitation choices."""
function phytoplankton_tendency(
    config::TendencyConfig{Growth,OrganicCycling,Zooplankton,Limitation}; plankton_idx::Int
) where {Growth,OrganicCycling,Zooplankton,Limitation}
    nutrients = config.nutrients

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        plankton = tracer_values.plankton
        resources = nutrient_resources(tracer_values, nutrients)
        half_saturation_arrays = half_saturation_parameters(parameters, nutrients)
        half_saturations = map(K -> K[plankton_idx], half_saturation_arrays)

        P = plankton(plankton_idx)
        growth = plankton_growth(
            Val(Growth),
            Val(Limitation),
            resources,
            P,
            tracer_values.PAR,
            parameters,
            half_saturations,
            plankton_idx,
        )
        grazing = grazing_loss_sum(parameters, plankton, P, plankton_idx)
        m = @inbounds parameters.linear_mortality[plankton_idx]
        mort = linear_loss(P, m)

        return growth - grazing - mort
    end

    return CompiledEquation(f)
end
