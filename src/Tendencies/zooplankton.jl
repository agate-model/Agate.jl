"""Zooplankton tendency with preferential grazing gain and mortality losses."""
function zooplankton_tendency(
    config::TendencyConfig{Formulation,:preferential_grazing,Limitation};
    plankton_idx::Int,
) where {Formulation,Limitation}
    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        plankton = tracer_values.plankton
        Z = plankton(plankton_idx)

        gain = grazing_gain_sum(parameters, plankton, Z, plankton_idx)

        m_lin = @inbounds parameters.linear_mortality[plankton_idx]
        m_quad = @inbounds parameters.quadratic_mortality[plankton_idx]

        return gain - linear_loss(Z, m_lin) - quadratic_loss(Z, m_quad)
    end

    return CompiledEquation(f)
end
