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

using ....Utils: sum_over, TendencyContext, tendency_views

using ...Sums:
    grazing_unassimilated_loss_sum,
    grazing_loss_sum,
    grazing_gain_sum,
    smith_uptake_sum,
    linear_mortality_sum,
    quadratic_mortality_sum

export nutrient_default, detritus_default, phytoplankton_default, zooplankton_default

@inline linear_mortality_loss(parameters, idx::Int, P) =
    linear_loss(P, parameters.linear_mortality[idx])

@inline quadratic_mortality_loss(parameters, idx::Int, P) =
    quadratic_loss(P, parameters.quadratic_mortality[idx])

@inline function smith_growth(parameters, idx::Int, N, P, PAR)
    return smith_single_nutrient_growth(
        N,
        P,
        PAR,
        parameters.maximum_growth_rate[idx],
        parameters.nutrient_half_saturation[idx],
        parameters.alpha[idx],
    )
end

"""Nutrient tendency with Smith growth and mortality/remineralization."""
function nutrient_default(plankton_syms)
    n_plankton = length(plankton_syms)

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
        tendency, parameters, vals = tendency_views(bgc, args)

        N = vals.N
        D = vals.D
        PAR = vals.PAR

        lin_sum = linear_mortality_sum(n_plankton, vals, parameters.linear_mortality)
        quad_sum = quadratic_mortality_sum(n_plankton, vals, parameters.quadratic_mortality)

        uptake = smith_uptake_sum(
            n_plankton,
            vals,
            N,
            PAR,
            parameters.maximum_growth_rate,
            parameters.nutrient_half_saturation,
            parameters.alpha,
        )

        export_frac = parameters.mortality_export_fraction
        remin_term = linear_remineralization(D, parameters.detritus_remineralization)

        return export_frac * (lin_sum + quad_sum) + remin_term - uptake
    end

    return CompiledEquation(f, requirements)
end

"""Detritus tendency from mortality, sloppy feeding, and remineralization."""
function detritus_default(plankton_syms)
    n_plankton = length(plankton_syms)

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
        tendency, parameters, vals = tendency_views(bgc, args)

        D = vals.D

        lin_sum = linear_mortality_sum(n_plankton, vals, parameters.linear_mortality)
        quad_sum = quadratic_mortality_sum(n_plankton, vals, parameters.quadratic_mortality)

        unassimilated = grazing_unassimilated_loss_sum(tendency)

        export_frac = parameters.mortality_export_fraction
        remin_term = linear_remineralization(D, parameters.detritus_remineralization)

        return (one(export_frac) - export_frac) * (lin_sum + quad_sum) + unassimilated - remin_term
    end

    return CompiledEquation(f, requirements)
end

"""Phytoplankton tendency with Smith growth, grazing loss, and linear mortality."""
function phytoplankton_default(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
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
        tendency, parameters, vals = tendency_views(bgc, args)

        N = vals.N
        P = vals.plankton(plankton_idx)
        PAR = vals.PAR

        growth = smith_growth(parameters, plankton_idx, N, P, PAR)

        grazing = grazing_loss_sum(tendency, P, plankton_idx)
        mort = linear_mortality_loss(parameters, plankton_idx, P)

        return growth - grazing - mort
    end

    return CompiledEquation(f, requirements)
end

"""Zooplankton tendency with preferential grazing gain and mortality losses."""
function zooplankton_default(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
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
        tendency, parameters, vals = tendency_views(bgc, args)

        Z = vals.plankton(plankton_idx)

        gain = grazing_gain_sum(tendency, Z, plankton_idx)

        lin = linear_mortality_loss(parameters, plankton_idx, Z)
        quad = quadratic_mortality_loss(parameters, plankton_idx, Z)

        return gain - lin - quad
    end

    return CompiledEquation(f, requirements)
end

end # module
