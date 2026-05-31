"""Organic matter tendency from plankton losses and linear remineralization."""
function organic_matter_tendency(
    config::TendencyConfig{:geider_dom_pom,Zooplankton,Limitation};
    target::Symbol,
    remineralization::Symbol,
    fraction::Symbol,
    stoichiometry::Symbol=:one,
) where {Zooplankton,Limitation}
    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        plankton = tracer_values.plankton
        pool = tracer_value(tracer_values, Val(target))

        mortality = mortality_loss_sum(
            plankton, parameters.linear_mortality, parameters.quadratic_mortality
        )
        grazing = grazing_unassimilated_loss_sum(parameters, plankton)
        remin = linear_remineralization(pool, parameter_value(parameters, Val(remineralization)))

        frac = organic_fraction(parameters, Val(fraction))
        ratio = stoichiometry_multiplier(parameters, Val(stoichiometry))

        return frac * ratio * (mortality + grazing) - remin
    end

    return CompiledEquation(f)
end
