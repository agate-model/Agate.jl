"""Simple detritus pool tendency for exported mortality and unassimilated grazing."""
function detritus_tendency(
    config::TendencyConfig{Growth,:simple_detritus,Zooplankton,Limitation};
    target::Symbol=:D,
    remineralization::Symbol=:detritus_remineralization,
    export_fraction::Symbol=:mortality_export_fraction,
) where {Growth,Zooplankton,Limitation}
    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        plankton = tracer_values.plankton
        D = tracer_value(tracer_values, Val(target))

        mortality = mortality_loss_sum(
            plankton, parameters.linear_mortality, parameters.quadratic_mortality
        )
        unassimilated = grazing_unassimilated_loss_sum(parameters, plankton)

        export_frac = parameter_value(parameters, Val(export_fraction))
        remin_term = linear_remineralization(D, parameter_value(parameters, Val(remineralization)))

        return (one(export_frac) - export_frac) * mortality + unassimilated - remin_term
    end

    return CompiledEquation(f)
end
