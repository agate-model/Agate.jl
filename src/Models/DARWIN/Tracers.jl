"""Tracer tendency functors for the simplified DARWIN model.

Predation terms use the canonical rectangular interaction matrices stored in
`bgc.parameters.interactions` (consumer-by-prey).

GPU notes
---------
The compiled tracer equations operate directly on positional tracer arguments
via `bgc.tracers`. No runtime Symbol indexing is performed in kernel-callable code.
"""

module Tracers

using ....Functors: CompiledEquation, req

using ....Library.Mortality: LinearLoss, QuadraticLoss
using ....Library.Photosynthesis: TwoNutrientGrowthGeider
using ....Library.Remineralization: LinearRemineralization

using ....Utils: sum_over, KernelBundle

using ...PredationSums: _grazing_assimilation_loss_sum, _grazing_loss_sum, _grazing_gain_sum

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

@inline function _uptake_sum_geider(
    p,
    tracers,
    args,
    npl::Int,
    DIN,
    PO4,
    PAR,
)
    sum_over(npl, zero(DIN)) do i
        P = tracers.plankton(args, i)
        TwoNutrientGrowthGeider(
            p.maximum_growth_rate[i],
            p.half_saturation_DIN[i],
            p.half_saturation_PO4[i],
            p.photosynthetic_slope[i],
            p.chlorophyll_to_carbon_ratio[i],
        )(
            DIN, PO4, P, PAR
        )
    end
end

@inline function _mortality_loss_sum(p, kb::KernelBundle, npl::Int)
    tracers = kb.tracers
    args = kb.args

    FT = eltype(p.maximum_growth_rate)
    z = zero(FT)

    sum_over(npl, z) do i
        P = tracers.plankton(args, i)
        LinearLoss(p.linear_mortality[i])(P) + QuadraticLoss(p.quadratic_mortality[i])(P)
    end
end

"""DIC tendency with Geider-style growth (carbon units)."""
function DIC_geider_light(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
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
        tracers = bgc.tracers

        DIN = tracers.DIN(args)
        PO4 = tracers.PO4(args)
        DOC = tracers.DOC(args)
        POC = tracers.POC(args)
        PAR = tracers.PAR(args)

        uptake = _uptake_sum_geider(p, tracers, args, npl, DIN, PO4, PAR)

        dic_remin =
            LinearRemineralization(p.DOC_remineralization)(DOC) +
            LinearRemineralization(p.POC_remineralization)(POC)

        return dic_remin - uptake
    end

    return CompiledEquation(f, requirements)
end

"""DIN tendency assuming fixed stoichiometry (N:C)."""
function DIN_geider_light(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
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
        tracers = bgc.tracers

        DIN = tracers.DIN(args)
        PO4 = tracers.PO4(args)
        DON = tracers.DON(args)
        PON = tracers.PON(args)
        PAR = tracers.PAR(args)

        uptake = _uptake_sum_geider(p, tracers, args, npl, DIN, PO4, PAR)

        din_remin =
            LinearRemineralization(p.DON_remineralization)(DON) +
            LinearRemineralization(p.PON_remineralization)(PON)

        return din_remin - p.nitrogen_to_carbon * uptake
    end

    return CompiledEquation(f, requirements)
end

"""PO4 tendency assuming fixed stoichiometry (P:C)."""
function PO4_geider_light(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
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
        tracers = bgc.tracers

        DIN = tracers.DIN(args)
        PO4 = tracers.PO4(args)
        DOP = tracers.DOP(args)
        POP = tracers.POP(args)
        PAR = tracers.PAR(args)

        uptake = _uptake_sum_geider(p, tracers, args, npl, DIN, PO4, PAR)

        po4_remin =
            LinearRemineralization(p.DOP_remineralization)(DOP) +
            LinearRemineralization(p.POP_remineralization)(POP)

        return po4_remin - p.phosphorus_to_carbon * uptake
    end

    return CompiledEquation(f, requirements)
end

# --- Organic matter ---------------------------------------------------------

"""DOC tendency from plankton losses and remineralization."""
function DOC_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
        scalars=(:DOM_POM_fractionation, :DOC_remineralization),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters
        tracers = bgc.tracers

        kb = KernelBundle(tracers, args)

        DOC = tracers.DOC(args)

        M = _mortality_loss_sum(p, kb, npl)
        g = _grazing_assimilation_loss_sum(p, kb)
        R = LinearRemineralization(p.DOC_remineralization)(DOC)

        frac = one(p.DOM_POM_fractionation) - p.DOM_POM_fractionation
        return frac * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""POC tendency from plankton losses and remineralization."""
function POC_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
        scalars=(:DOM_POM_fractionation, :POC_remineralization),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters
        tracers = bgc.tracers

        kb = KernelBundle(tracers, args)

        POC = tracers.POC(args)

        M = _mortality_loss_sum(p, kb, npl)
        g = _grazing_assimilation_loss_sum(p, kb)
        R = LinearRemineralization(p.POC_remineralization)(POC)

        return p.DOM_POM_fractionation * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""DON tendency assuming fixed stoichiometry (N:C)."""
function DON_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
        scalars=(:DOM_POM_fractionation, :DON_remineralization, :nitrogen_to_carbon),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters
        tracers = bgc.tracers

        kb = KernelBundle(tracers, args)

        DON = tracers.DON(args)

        M = _mortality_loss_sum(p, kb, npl)
        g = _grazing_assimilation_loss_sum(p, kb)
        R = LinearRemineralization(p.DON_remineralization)(DON)

        frac = one(p.DOM_POM_fractionation) - p.DOM_POM_fractionation
        return frac * p.nitrogen_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""PON tendency assuming fixed stoichiometry (N:C)."""
function PON_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
        scalars=(:DOM_POM_fractionation, :PON_remineralization, :nitrogen_to_carbon),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters
        tracers = bgc.tracers

        kb = KernelBundle(tracers, args)

        PON = tracers.PON(args)

        M = _mortality_loss_sum(p, kb, npl)
        g = _grazing_assimilation_loss_sum(p, kb)
        R = LinearRemineralization(p.PON_remineralization)(PON)

        return p.DOM_POM_fractionation * p.nitrogen_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""DOP tendency assuming fixed stoichiometry (P:C)."""
function DOP_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
        scalars=(:DOM_POM_fractionation, :DOP_remineralization, :phosphorus_to_carbon),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters
        tracers = bgc.tracers

        kb = KernelBundle(tracers, args)

        DOP = tracers.DOP(args)

        M = _mortality_loss_sum(p, kb, npl)
        g = _grazing_assimilation_loss_sum(p, kb)
        R = LinearRemineralization(p.DOP_remineralization)(DOP)

        frac = one(p.DOM_POM_fractionation) - p.DOM_POM_fractionation
        return frac * p.phosphorus_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""POP tendency assuming fixed stoichiometry (P:C)."""
function POP_default(plankton_syms)
    npl = length(plankton_syms)

    requirements = req(;
        scalars=(:DOM_POM_fractionation, :POP_remineralization, :phosphorus_to_carbon),
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters
        tracers = bgc.tracers

        kb = KernelBundle(tracers, args)

        POP = tracers.POP(args)

        M = _mortality_loss_sum(p, kb, npl)
        g = _grazing_assimilation_loss_sum(p, kb)
        R = LinearRemineralization(p.POP_remineralization)(POP)

        return p.DOM_POM_fractionation * p.phosphorus_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

# --- Plankton ---------------------------------------------------------------

"""Phytoplankton tendency with Geider-style, two-nutrient growth."""
function phytoplankton_growth_two_nutrients_geider_light(
    plankton_syms,
    plankton_sym::Symbol,
    plankton_idx::Int,
)
    requirements = req(;
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
        tracers = bgc.tracers

        kb = KernelBundle(tracers, args)

        DIN = tracers.DIN(args)
        PO4 = tracers.PO4(args)
        P = tracers.plankton(args, plankton_idx)
        PAR = tracers.PAR(args)

        growth = TwoNutrientGrowthGeider(
            p.maximum_growth_rate[plankton_idx],
            p.half_saturation_DIN[plankton_idx],
            p.half_saturation_PO4[plankton_idx],
            p.photosynthetic_slope[plankton_idx],
            p.chlorophyll_to_carbon_ratio[plankton_idx],
        )(
            DIN, PO4, P, PAR
        )

        grazing = _grazing_loss_sum(p, kb, P, plankton_idx)

        mort = LinearLoss(p.linear_mortality[plankton_idx])(P)

        return growth - grazing - mort
    end

    return CompiledEquation(f, requirements)
end

"""Zooplankton tendency with preferential grazing gain."""
function zooplankton_default(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    requirements = req(;
        vectors=(
            :linear_mortality,
            :quadratic_mortality,
            :maximum_predation_rate,
            :holling_half_saturation,
        ),
        matrices=(:palatability_matrix, :assimilation_matrix),
    )

    f = function (bgc, x, y, z, t, args...)
        p = bgc.parameters
        tracers = bgc.tracers

        kb = KernelBundle(tracers, args)

        Z = tracers.plankton(args, plankton_idx)

        gain = _grazing_gain_sum(p, kb, Z, plankton_idx)

        lin = LinearLoss(p.linear_mortality[plankton_idx])(Z)
        quad = QuadraticLoss(p.quadratic_mortality[plankton_idx])(Z)

        return gain - lin - quad
    end

    return CompiledEquation(f, requirements)
end

end # module
