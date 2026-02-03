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

using ....Library.Mortality: LinearLoss, QuadraticLoss
using ....Library.Photosynthesis: SingleNutrientGrowthSmith
using ....Library.Remineralization: LinearRemineralization

using ....Utils: sum_over, TendencyContext

using ...PredationSums: _grazing_assimilation_loss_sum, _grazing_loss_sum, _grazing_gain_sum

export nutrient_default, detritus_default, phytoplankton_default, zooplankton_default

@inline function _uptake_sum_smith(tendency::TendencyContext, n_plankton::Int, N, PAR)
    parameters = tendency.parameters
    tracers = tendency.tracers
    args = tendency.args

    sum_over(n_plankton, zero(N)) do i
        P = tracers.plankton(args, i)
        SingleNutrientGrowthSmith(
            parameters.maximum_growth_rate[i],
            parameters.nutrient_half_saturation[i],
            parameters.alpha[i],
        )(N, P, PAR)
    end
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        N = tracers.N(args)
        D = tracers.D(args)
        PAR = tracers.PAR(args)

        lin_sum = sum_over(n_plankton, zero(N)) do i
            P = tracers.plankton(args, i)
            LinearLoss(parameters.linear_mortality[i])(P)
        end

        quad_sum = sum_over(n_plankton, zero(N)) do i
            P = tracers.plankton(args, i)
            QuadraticLoss(parameters.quadratic_mortality[i])(P)
        end

        uptake = _uptake_sum_smith(tendency, n_plankton, N, PAR)

        export_frac = parameters.mortality_export_fraction
        remin = LinearRemineralization(parameters.detritus_remineralization)(D)

        return export_frac * (lin_sum + quad_sum) + remin - uptake
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        D = tracers.D(args)

        lin_sum = sum_over(n_plankton, zero(D)) do i
            P = tracers.plankton(args, i)
            LinearLoss(parameters.linear_mortality[i])(P)
        end

        quad_sum = sum_over(n_plankton, zero(D)) do i
            P = tracers.plankton(args, i)
            QuadraticLoss(parameters.quadratic_mortality[i])(P)
        end

        assim_loss = _grazing_assimilation_loss_sum(tendency)

        export_frac = parameters.mortality_export_fraction
        remin = LinearRemineralization(parameters.detritus_remineralization)(D)

        return (one(export_frac) - export_frac) * (lin_sum + quad_sum) + assim_loss - remin
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        N = tracers.N(args)
        P = tracers.plankton(args, plankton_idx)
        PAR = tracers.PAR(args)

        growth = SingleNutrientGrowthSmith(
            parameters.maximum_growth_rate[plankton_idx],
            parameters.nutrient_half_saturation[plankton_idx],
            parameters.alpha[plankton_idx],
        )(N, P, PAR)

        grazing = _grazing_loss_sum(tendency, P, plankton_idx)
        mort = LinearLoss(parameters.linear_mortality[plankton_idx])(P)

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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        Z = tracers.plankton(args, plankton_idx)

        gain = _grazing_gain_sum(tendency, Z, plankton_idx)

        lin = LinearLoss(parameters.linear_mortality[plankton_idx])(Z)
        quad = QuadraticLoss(parameters.quadratic_mortality[plankton_idx])(Z)

        return gain - lin - quad
    end

    return CompiledEquation(f, requirements)
end

end # module
