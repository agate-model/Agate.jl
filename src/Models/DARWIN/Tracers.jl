"""Tracer tendency functors for the simplified DARWIN model."""

module Tracers

using ....Functors: CompiledEquation, req

using ....Library.Mortality: LinearLoss, QuadraticLoss
using ....Library.Photosynthesis: TwoNutrientGrowthGeider
using ....Library.Predation: PreferentialPredationLoss, PreferentialPredationGain, PreferentialPredationAssimilationLoss
using ....Library.Remineralization: LinearRemineralization

export DIC_geider_light,
    DIN_geider_light,
    PO4_geider_light,
    POC_default,
    DOC_default,
    PON_default,
    DON_default,
    POP_default,
    DOP_default,
    phytoplankton_growth_two_nutrients_geider_light,
    zooplankton_default


const _N_BIO_TRACERS = 9

"""DIC tendency with Geider-style growth (carbon units)."""
function DIC_geider_light(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOC_remineralization, :POC_remineralization),
        vectors=(
            :maximum_growth_rate,
            :half_saturation_DIN,
            :half_saturation_PO4,
            :photosynthetic_slope,
            :chlorophyll_to_carbon_ratio,
        ),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        DIN = args[2]
        PO4 = args[3]
        DOC = args[4]
        POC = args[5]
        PAR = args[_N_BIO_TRACERS + npl + 1]

        growth_sum = zero(DIN)
        @inbounds for i in 1:npl
            P = args[_N_BIO_TRACERS + i]
            growth_sum += TwoNutrientGrowthGeider(
                p.maximum_growth_rate[i],
                p.half_saturation_DIN[i],
                p.half_saturation_PO4[i],
                p.photosynthetic_slope[i],
                p.chlorophyll_to_carbon_ratio[i],
            )(DIN, PO4, P, PAR)
        end

        dic_remin = LinearRemineralization(p.DOC_remineralization)(DOC) +
                    LinearRemineralization(p.POC_remineralization)(POC)

        return dic_remin - growth_sum
    end

    return CompiledEquation(f, requirements)
end

"""DIN tendency assuming fixed stoichiometry (N:C)."""
function DIN_geider_light(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DON_remineralization, :PON_remineralization, :nitrogen_to_carbon),
        vectors=(
            :maximum_growth_rate,
            :half_saturation_DIN,
            :half_saturation_PO4,
            :photosynthetic_slope,
            :chlorophyll_to_carbon_ratio,
        ),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        DIN = args[2]
        PO4 = args[3]
        DON = args[6]
        PON = args[7]
        PAR = args[_N_BIO_TRACERS + npl + 1]

        growth_sum = zero(DIN)
        @inbounds for i in 1:npl
            P = args[_N_BIO_TRACERS + i]
            growth_sum += TwoNutrientGrowthGeider(
                p.maximum_growth_rate[i],
                p.half_saturation_DIN[i],
                p.half_saturation_PO4[i],
                p.photosynthetic_slope[i],
                p.chlorophyll_to_carbon_ratio[i],
            )(DIN, PO4, P, PAR)
        end

        din_remin = LinearRemineralization(p.DON_remineralization)(DON) +
                    LinearRemineralization(p.PON_remineralization)(PON)

        return din_remin - p.nitrogen_to_carbon * growth_sum
    end

    return CompiledEquation(f, requirements)
end

"""PO4 tendency assuming fixed stoichiometry (P:C)."""
function PO4_geider_light(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOP_remineralization, :POP_remineralization, :phosphorus_to_carbon),
        vectors=(
            :maximum_growth_rate,
            :half_saturation_DIN,
            :half_saturation_PO4,
            :photosynthetic_slope,
            :chlorophyll_to_carbon_ratio,
        ),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        DIN = args[2]
        PO4 = args[3]
        DOP = args[8]
        POP = args[9]
        PAR = args[_N_BIO_TRACERS + npl + 1]

        growth_sum = zero(DIN)
        @inbounds for i in 1:npl
            P = args[_N_BIO_TRACERS + i]
            growth_sum += TwoNutrientGrowthGeider(
                p.maximum_growth_rate[i],
                p.half_saturation_DIN[i],
                p.half_saturation_PO4[i],
                p.photosynthetic_slope[i],
                p.chlorophyll_to_carbon_ratio[i],
            )(DIN, PO4, P, PAR)
        end

        po4_remin = LinearRemineralization(p.DOP_remineralization)(DOP) +
                    LinearRemineralization(p.POP_remineralization)(POP)

        return po4_remin - p.phosphorus_to_carbon * growth_sum
    end

    return CompiledEquation(f, requirements)
end

# --- Organic matter ---------------------------------------------------------

"""DOC tendency from plankton losses and remineralization."""
function DOC_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :DOC_remineralization),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        DOC = args[4]

        lin_sum = zero(DOC)
        quad_sum = zero(DOC)
        @inbounds for i in 1:npl
            P = args[_N_BIO_TRACERS + i]
            lin_sum += LinearLoss(p.linear_mortality[i])(P)
            quad_sum += QuadraticLoss(p.quadratic_mortality[i])(P)
        end

        assim_loss = zero(DOC)
        @inbounds for predator_idx in 1:npl
            predator = args[_N_BIO_TRACERS + predator_idx]
            gmax = p.maximum_predation_rate[predator_idx]
            K = p.holling_half_saturation[predator_idx]
            for prey_idx in 1:npl
                prey = args[_N_BIO_TRACERS + prey_idx]
                β = p.assimilation_matrix[predator_idx, prey_idx]
                ϕ = p.palatability_matrix[predator_idx, prey_idx]
                assim_loss += PreferentialPredationAssimilationLoss(β, gmax, K, ϕ)(prey, predator)
            end
        end

        base = lin_sum + quad_sum + assim_loss
        remin = LinearRemineralization(p.DOC_remineralization)(DOC)

        return (one(p.DOM_POM_fractionation) - p.DOM_POM_fractionation) * base - remin
    end

    return CompiledEquation(f, requirements)
end

"""POC tendency from plankton losses and remineralization."""
function POC_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :POC_remineralization),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        POC = args[5]

        lin_sum = zero(POC)
        quad_sum = zero(POC)
        @inbounds for i in 1:npl
            P = args[_N_BIO_TRACERS + i]
            lin_sum += LinearLoss(p.linear_mortality[i])(P)
            quad_sum += QuadraticLoss(p.quadratic_mortality[i])(P)
        end

        assim_loss = zero(POC)
        @inbounds for predator_idx in 1:npl
            predator = args[_N_BIO_TRACERS + predator_idx]
            gmax = p.maximum_predation_rate[predator_idx]
            K = p.holling_half_saturation[predator_idx]
            for prey_idx in 1:npl
                prey = args[_N_BIO_TRACERS + prey_idx]
                β = p.assimilation_matrix[predator_idx, prey_idx]
                ϕ = p.palatability_matrix[predator_idx, prey_idx]
                assim_loss += PreferentialPredationAssimilationLoss(β, gmax, K, ϕ)(prey, predator)
            end
        end

        base = lin_sum + quad_sum + assim_loss
        remin = LinearRemineralization(p.POC_remineralization)(POC)

        return p.DOM_POM_fractionation * base - remin
    end

    return CompiledEquation(f, requirements)
end

"""DON tendency assuming fixed stoichiometry (N:C)."""
function DON_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :DON_remineralization, :nitrogen_to_carbon),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        DON = args[6]

        lin_sum = zero(DON)
        quad_sum = zero(DON)
        @inbounds for i in 1:npl
            P = args[_N_BIO_TRACERS + i]
            lin_sum += LinearLoss(p.linear_mortality[i])(P)
            quad_sum += QuadraticLoss(p.quadratic_mortality[i])(P)
        end

        assim_loss = zero(DON)
        @inbounds for predator_idx in 1:npl
            predator = args[_N_BIO_TRACERS + predator_idx]
            gmax = p.maximum_predation_rate[predator_idx]
            K = p.holling_half_saturation[predator_idx]
            for prey_idx in 1:npl
                prey = args[_N_BIO_TRACERS + prey_idx]
                β = p.assimilation_matrix[predator_idx, prey_idx]
                ϕ = p.palatability_matrix[predator_idx, prey_idx]
                assim_loss += PreferentialPredationAssimilationLoss(β, gmax, K, ϕ)(prey, predator)
            end
        end

        base = lin_sum + quad_sum + assim_loss
        remin = LinearRemineralization(p.DON_remineralization)(DON)

        frac = one(p.DOM_POM_fractionation) - p.DOM_POM_fractionation
        return frac * p.nitrogen_to_carbon * base - remin
    end

    return CompiledEquation(f, requirements)
end

"""PON tendency assuming fixed stoichiometry (N:C)."""
function PON_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :PON_remineralization, :nitrogen_to_carbon),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        PON = args[7]

        lin_sum = zero(PON)
        quad_sum = zero(PON)
        @inbounds for i in 1:npl
            P = args[_N_BIO_TRACERS + i]
            lin_sum += LinearLoss(p.linear_mortality[i])(P)
            quad_sum += QuadraticLoss(p.quadratic_mortality[i])(P)
        end

        assim_loss = zero(PON)
        @inbounds for predator_idx in 1:npl
            predator = args[_N_BIO_TRACERS + predator_idx]
            gmax = p.maximum_predation_rate[predator_idx]
            K = p.holling_half_saturation[predator_idx]
            for prey_idx in 1:npl
                prey = args[_N_BIO_TRACERS + prey_idx]
                β = p.assimilation_matrix[predator_idx, prey_idx]
                ϕ = p.palatability_matrix[predator_idx, prey_idx]
                assim_loss += PreferentialPredationAssimilationLoss(β, gmax, K, ϕ)(prey, predator)
            end
        end

        base = lin_sum + quad_sum + assim_loss
        remin = LinearRemineralization(p.PON_remineralization)(PON)

        return p.DOM_POM_fractionation * p.nitrogen_to_carbon * base - remin
    end

    return CompiledEquation(f, requirements)
end

"""DOP tendency assuming fixed stoichiometry (P:C)."""
function DOP_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :DOP_remineralization, :phosphorus_to_carbon),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        DOP = args[8]

        lin_sum = zero(DOP)
        quad_sum = zero(DOP)
        @inbounds for i in 1:npl
            P = args[_N_BIO_TRACERS + i]
            lin_sum += LinearLoss(p.linear_mortality[i])(P)
            quad_sum += QuadraticLoss(p.quadratic_mortality[i])(P)
        end

        assim_loss = zero(DOP)
        @inbounds for predator_idx in 1:npl
            predator = args[_N_BIO_TRACERS + predator_idx]
            gmax = p.maximum_predation_rate[predator_idx]
            K = p.holling_half_saturation[predator_idx]
            for prey_idx in 1:npl
                prey = args[_N_BIO_TRACERS + prey_idx]
                β = p.assimilation_matrix[predator_idx, prey_idx]
                ϕ = p.palatability_matrix[predator_idx, prey_idx]
                assim_loss += PreferentialPredationAssimilationLoss(β, gmax, K, ϕ)(prey, predator)
            end
        end

        base = lin_sum + quad_sum + assim_loss
        remin = LinearRemineralization(p.DOP_remineralization)(DOP)

        frac = one(p.DOM_POM_fractionation) - p.DOM_POM_fractionation
        return frac * p.phosphorus_to_carbon * base - remin
    end

    return CompiledEquation(f, requirements)
end

"""POP tendency assuming fixed stoichiometry (P:C)."""
function POP_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(
        scalars=(:DOM_POM_fractionation, :POP_remineralization, :phosphorus_to_carbon),
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        POP = args[9]

        lin_sum = zero(POP)
        quad_sum = zero(POP)
        @inbounds for i in 1:npl
            P = args[_N_BIO_TRACERS + i]
            lin_sum += LinearLoss(p.linear_mortality[i])(P)
            quad_sum += QuadraticLoss(p.quadratic_mortality[i])(P)
        end

        assim_loss = zero(POP)
        @inbounds for predator_idx in 1:npl
            predator = args[_N_BIO_TRACERS + predator_idx]
            gmax = p.maximum_predation_rate[predator_idx]
            K = p.holling_half_saturation[predator_idx]
            for prey_idx in 1:npl
                prey = args[_N_BIO_TRACERS + prey_idx]
                β = p.assimilation_matrix[predator_idx, prey_idx]
                ϕ = p.palatability_matrix[predator_idx, prey_idx]
                assim_loss += PreferentialPredationAssimilationLoss(β, gmax, K, ϕ)(prey, predator)
            end
        end

        base = lin_sum + quad_sum + assim_loss
        remin = LinearRemineralization(p.POP_remineralization)(POP)

        return p.DOM_POM_fractionation * p.phosphorus_to_carbon * base - remin
    end

    return CompiledEquation(f, requirements)
end

# --- Plankton ---------------------------------------------------------------

"""Phytoplankton tendency with Geider-style, two-nutrient growth."""
function phytoplankton_growth_two_nutrients_geider_light(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    npl = length(plankton_syms)

    requirements = req(
        vectors=(
            :maximum_growth_rate,
            :half_saturation_DIN,
            :half_saturation_PO4,
            :photosynthetic_slope,
            :chlorophyll_to_carbon_ratio,
            :linear_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix,),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        DIN = args[2]
        PO4 = args[3]
        P = args[_N_BIO_TRACERS + plankton_idx]
        PAR = args[_N_BIO_TRACERS + npl + 1]

        growth = TwoNutrientGrowthGeider(
            p.maximum_growth_rate[plankton_idx],
            p.half_saturation_DIN[plankton_idx],
            p.half_saturation_PO4[plankton_idx],
            p.photosynthetic_slope[plankton_idx],
            p.chlorophyll_to_carbon_ratio[plankton_idx],
        )(DIN, PO4, P, PAR)

        grazing = zero(P)
        @inbounds for predator_idx in 1:npl
            predator = args[_N_BIO_TRACERS + predator_idx]
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

"""Zooplankton tendency with preferential grazing gain."""
function zooplankton_default(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    npl = length(plankton_syms)

    requirements = req(
        vectors=(:linear_mortality, :quadratic_mortality, :maximum_predation_rate, :holling_half_saturation),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters

        Z = args[_N_BIO_TRACERS + plankton_idx]

        gmax = p.maximum_predation_rate[plankton_idx]
        K = p.holling_half_saturation[plankton_idx]

        gain = zero(Z)
        @inbounds for prey_idx in 1:npl
            prey = args[_N_BIO_TRACERS + prey_idx]
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
