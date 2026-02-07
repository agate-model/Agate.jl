"""Tracer tendency functors for the NiPiZD model.

Predation terms use the canonical rectangular interaction matrices stored in
`bgc.parameters.interactions` (consumer-by-prey).

GPU notes
---------
The compiled tracer equations operate directly on positional tracer arguments
via `bgc.tracers`. No runtime Symbol indexing is performed in kernel-callable
code.
"""

module Tracers

using ....Functors: CompiledEquation, Requirements

using ....Library.Mortality: linear_loss, quadratic_loss
using ....Library.Photosynthesis: smith_single_nutrient_growth
using ....Library.Remineralization: linear_remineralization

using ....Utils: tendency_inputs

using ...Sums:
    grazing_unassimilated_loss_sum,
    grazing_loss_sum,
    grazing_gain_sum,
    smith_uptake_sum,
    mortality_loss_sum

export nutrient_default, detritus_default, phytoplankton_default, zooplankton_default

"""Nutrient tendency with Smith growth and mortality/remineralization."""
function nutrient_default()
    requirements = Requirements(;
        scalars=(:detritus_remineralization, :mortality_export_fraction),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_growth_rate,
            :nutrient_half_saturation,
            :alpha,
        ),
    )

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        plankton = tracer_values.plankton

        N = tracer_values.N
        D = tracer_values.D
        PAR = tracer_values.PAR

        mortality = mortality_loss_sum(
            plankton, parameters.linear_mortality, parameters.quadratic_mortality
        )

        uptake = smith_uptake_sum(
            plankton,
            N,
            PAR,
            parameters.maximum_growth_rate,
            parameters.nutrient_half_saturation,
            parameters.alpha,
        )

        export_frac = parameters.mortality_export_fraction
        remin_term = linear_remineralization(D, parameters.detritus_remineralization)

        return export_frac * mortality + remin_term - uptake
    end

    return CompiledEquation(f, requirements)
end

"""Detritus tendency from mortality, sloppy feeding, and remineralization."""
function detritus_default()
    requirements = Requirements(;
        scalars=(:detritus_remineralization, :mortality_export_fraction),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        plankton = tracer_values.plankton
        D = tracer_values.D

        mortality = mortality_loss_sum(
            plankton, parameters.linear_mortality, parameters.quadratic_mortality
        )

        unassimilated = grazing_unassimilated_loss_sum(parameters, plankton)

        export_frac = parameters.mortality_export_fraction
        remin_term = linear_remineralization(D, parameters.detritus_remineralization)

        return (one(export_frac) - export_frac) * mortality + unassimilated - remin_term
    end

    return CompiledEquation(f, requirements)
end

"""Phytoplankton tendency with Smith growth, grazing loss, and linear mortality."""
function phytoplankton_default(plankton_idx::Int)
    requirements = Requirements(;
        vectors=(
            :maximum_growth_rate,
            :nutrient_half_saturation,
            :alpha,
            :linear_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix,),
    )

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        plankton = tracer_values.plankton

        N = tracer_values.N
        PAR = tracer_values.PAR
        P = plankton(plankton_idx)

        μmax = @inbounds parameters.maximum_growth_rate[plankton_idx]
        K = @inbounds parameters.nutrient_half_saturation[plankton_idx]
        α = @inbounds parameters.alpha[plankton_idx]

        growth = smith_single_nutrient_growth(N, P, PAR, μmax, K, α)

        grazing = grazing_loss_sum(parameters, plankton, P, plankton_idx)

        m = @inbounds parameters.linear_mortality[plankton_idx]
        mort = linear_loss(P, m)

        return growth - grazing - mort
    end

    return CompiledEquation(f, requirements)
end

"""Zooplankton tendency with preferential grazing gain and mortality losses."""
function zooplankton_default(plankton_idx::Int)
    requirements = Requirements(;
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        parameters, tracer_values = tendency_inputs(bgc, args)

        plankton = tracer_values.plankton
        Z = plankton(plankton_idx)

        gain = grazing_gain_sum(parameters, plankton, Z, plankton_idx)

        m_lin = @inbounds parameters.linear_mortality[plankton_idx]
        m_quad = @inbounds parameters.quadratic_mortality[plankton_idx]

        lin = linear_loss(Z, m_lin)
        quad = quadratic_loss(Z, m_quad)

        return gain - lin - quad
    end

    return CompiledEquation(f, requirements)
end

end # module
