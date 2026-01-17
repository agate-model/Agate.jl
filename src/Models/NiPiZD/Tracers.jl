"""Tracer tendency functors for the NiPiZD model."""

module Tracers

using ....Equations: req
using ....Functors: CompiledEquation

using ....Library.Mortality: LinearLoss, QuadraticLoss
using ....Library.Photosynthesis: SingleNutrientGrowthSmith
using ....Library.Predation: PreferentialPredationLoss, PreferentialPredationGain, PreferentialPredationAssimilationLoss
using ....Library.Remineralization: LinearRemineralization

export nutrient_default, detritus_default, phytoplankton_default, zooplankton_default

@inline _pview(bgc) = hasproperty(bgc.parameters, :data) ? getproperty(bgc.parameters, :data) : bgc.parameters

"""Nutrient tendency with Smith growth and mortality/remineralization."""
function nutrient_default(PV, plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:detritus_remineralization, :mortality_export_fraction),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_growth_rate, :nutrient_half_saturation, :alpha),
    )

    f = function (bgc, x, y, z, t, args...)
        p = _pview(bgc)

        N = args[1]
        D = args[2]
        PAR = args[2 + npl + 1]

        lin_sum = zero(N)
        quad_sum = zero(N)
        growth_sum = zero(N)

        @inbounds for i in 1:npl
            P = args[2 + i]
            lin_sum += LinearLoss(p.linear_mortality[i])(P)
            quad_sum += QuadraticLoss(p.quadratic_mortality[i])(P)
            growth_sum += SingleNutrientGrowthSmith(
                p.maximum_growth_rate[i],
                p.nutrient_half_saturation[i],
                p.alpha[i],
            )(N, P, PAR)
        end

        export_frac = p.mortality_export_fraction
        remin = LinearRemineralization(p.detritus_remineralization)(D)

        return export_frac * (lin_sum + quad_sum) + remin - growth_sum
    end

    return CompiledEquation(f, requirements)
end

"""Detritus tendency from mortality, sloppy feeding, and remineralization."""
function detritus_default(PV, plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
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
        p = _pview(bgc)

        D = args[2]

        lin_sum = zero(D)
        quad_sum = zero(D)

        @inbounds for i in 1:npl
            P = args[2 + i]
            lin_sum += LinearLoss(p.linear_mortality[i])(P)
            quad_sum += QuadraticLoss(p.quadratic_mortality[i])(P)
        end

        assim_loss = zero(D)
        @inbounds for predator_idx in 1:npl
            predator = args[2 + predator_idx]
            gmax = p.maximum_predation_rate[predator_idx]
            K = p.holling_half_saturation[predator_idx]
            for prey_idx in 1:npl
                prey = args[2 + prey_idx]
                β = p.assimilation_matrix[predator_idx, prey_idx]
                ϕ = p.palatability_matrix[predator_idx, prey_idx]
                assim_loss += PreferentialPredationAssimilationLoss(β, gmax, K, ϕ)(prey, predator)
            end
        end

        export_frac = p.mortality_export_fraction
        remin = LinearRemineralization(p.detritus_remineralization)(D)

        return (one(export_frac) - export_frac) * (lin_sum + quad_sum) + assim_loss - remin
    end

    return CompiledEquation(f, requirements)
end

"""Phytoplankton tendency with Smith growth, grazing loss, and linear mortality."""
function phytoplankton_default(PV, plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    npl = length(plankton_syms)

    requirements = req(
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
        p = _pview(bgc)

        N = args[1]
        P = args[2 + plankton_idx]
        PAR = args[2 + npl + 1]

        growth = SingleNutrientGrowthSmith(
            p.maximum_growth_rate[plankton_idx],
            p.nutrient_half_saturation[plankton_idx],
            p.alpha[plankton_idx],
        )(N, P, PAR)

        grazing = zero(P)
        @inbounds for predator_idx in 1:npl
            predator = args[2 + predator_idx]
            gmax = p.maximum_predation_rate[predator_idx]
            K = p.holling_half_saturation[predator_idx]
            ϕ = p.palatability_matrix[predator_idx, plankton_idx]
            grazing += PreferentialPredationLoss(gmax, K, ϕ)(P, predator)
        end

        mort = LinearLoss(p.linear_mortality[plankton_idx])(P)

        return growth - grazing - mort
    end

    return CompiledEquation(f, requirements)
end

"""Zooplankton tendency with preferential grazing gain and mortality losses."""
function zooplankton_default(PV, plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    npl = length(plankton_syms)

    requirements = req(
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = _pview(bgc)

        Z = args[2 + plankton_idx]

        gmax = p.maximum_predation_rate[plankton_idx]
        K = p.holling_half_saturation[plankton_idx]

        gain = zero(Z)
        @inbounds for prey_idx in 1:npl
            prey = args[2 + prey_idx]
            β = p.assimilation_matrix[plankton_idx, prey_idx]
            ϕ = p.palatability_matrix[plankton_idx, prey_idx]
            gain += PreferentialPredationGain(β, gmax, K, ϕ)(prey, Z)
        end

        lin = LinearLoss(p.linear_mortality[plankton_idx])(Z)
        quad = QuadraticLoss(p.quadratic_mortality[plankton_idx])(Z)

        return gain - lin - quad
    end

    return CompiledEquation(f, requirements)
end

end # module
