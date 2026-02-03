"""Tracer tendency functors for the simplified DARWIN model.

Predation terms use the canonical rectangular interaction matrices stored in
`bgc.parameters.interactions` (consumer-by-prey).

GPU notes
---------
The compiled tracer equations operate directly on positional tracer arguments
via `bgc.tracers`. No runtime Symbol indexing is performed in kernel-callable code.
"""

module Tracers

using ....Functors: CompiledEquation, Requirements

using ....Library.Mortality: LinearLoss, QuadraticLoss
using ....Library.Photosynthesis: TwoNutrientGrowthGeider
using ....Library.Remineralization: LinearRemineralization

using ....Utils: sum_over, TendencyContext

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

@inline function _uptake_sum_geider(tendency::TendencyContext, n_plankton::Int, DIN, PO4, PAR)
    parameters = tendency.parameters
    tracers = tendency.tracers
    args = tendency.args

    sum_over(n_plankton, zero(DIN)) do i
        P = tracers.plankton(args, i)
        TwoNutrientGrowthGeider(
            parameters.maximum_growth_rate[i],
            parameters.half_saturation_DIN[i],
            parameters.half_saturation_PO4[i],
            parameters.photosynthetic_slope[i],
            parameters.chlorophyll_to_carbon_ratio[i],
        )(
            DIN, PO4, P, PAR
        )
    end
end

@inline function _mortality_loss_sum(tendency::TendencyContext, n_plankton::Int)
    parameters = tendency.parameters
    tracers = tendency.tracers
    args = tendency.args

    FT = eltype(parameters.maximum_growth_rate)
    z = zero(FT)

    sum_over(n_plankton, z) do i
        P = tracers.plankton(args, i)
        LinearLoss(parameters.linear_mortality[i])(P) + QuadraticLoss(parameters.quadratic_mortality[i])(P)
    end
end

"""DIC tendency with Geider-style growth (carbon units)."""
function DIC_geider_light(plankton_syms)
    n_plankton = length(plankton_syms)

    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        DIN = tracers.DIN(args)
        PO4 = tracers.PO4(args)
        DOC = tracers.DOC(args)
        POC = tracers.POC(args)
        PAR = tracers.PAR(args)

        uptake = _uptake_sum_geider(tendency, n_plankton, DIN, PO4, PAR)

        dic_remin =
            LinearRemineralization(parameters.DOC_remineralization)(DOC) +
            LinearRemineralization(parameters.POC_remineralization)(POC)

        return dic_remin - uptake
    end

    return CompiledEquation(f, requirements)
end

"""DIN tendency assuming fixed stoichiometry (N:C)."""
function DIN_geider_light(plankton_syms)
    n_plankton = length(plankton_syms)

    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        DIN = tracers.DIN(args)
        PO4 = tracers.PO4(args)
        DON = tracers.DON(args)
        PON = tracers.PON(args)
        PAR = tracers.PAR(args)

        uptake = _uptake_sum_geider(tendency, n_plankton, DIN, PO4, PAR)

        din_remin =
            LinearRemineralization(parameters.DON_remineralization)(DON) +
            LinearRemineralization(parameters.PON_remineralization)(PON)

        return din_remin - parameters.nitrogen_to_carbon * uptake
    end

    return CompiledEquation(f, requirements)
end

"""PO4 tendency assuming fixed stoichiometry (P:C)."""
function PO4_geider_light(plankton_syms)
    n_plankton = length(plankton_syms)

    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        DIN = tracers.DIN(args)
        PO4 = tracers.PO4(args)
        DOP = tracers.DOP(args)
        POP = tracers.POP(args)
        PAR = tracers.PAR(args)

        uptake = _uptake_sum_geider(tendency, n_plankton, DIN, PO4, PAR)

        po4_remin =
            LinearRemineralization(parameters.DOP_remineralization)(DOP) +
            LinearRemineralization(parameters.POP_remineralization)(POP)

        return po4_remin - parameters.phosphorus_to_carbon * uptake
    end

    return CompiledEquation(f, requirements)
end

# --- Organic matter ---------------------------------------------------------

"""DOC tendency from plankton losses and remineralization."""
function DOC_default(plankton_syms)
    n_plankton = length(plankton_syms)

    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        DOC = tracers.DOC(args)

        M = _mortality_loss_sum(tendency, n_plankton)
        g = _grazing_assimilation_loss_sum(tendency)
        R = LinearRemineralization(parameters.DOC_remineralization)(DOC)

        frac = one(parameters.DOM_POM_fractionation) - parameters.DOM_POM_fractionation
        return frac * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""POC tendency from plankton losses and remineralization."""
function POC_default(plankton_syms)
    n_plankton = length(plankton_syms)

    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        POC = tracers.POC(args)

        M = _mortality_loss_sum(tendency, n_plankton)
        g = _grazing_assimilation_loss_sum(tendency)
        R = LinearRemineralization(parameters.POC_remineralization)(POC)

        return parameters.DOM_POM_fractionation * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""DON tendency assuming fixed stoichiometry (N:C)."""
function DON_default(plankton_syms)
    n_plankton = length(plankton_syms)

    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        DON = tracers.DON(args)

        M = _mortality_loss_sum(tendency, n_plankton)
        g = _grazing_assimilation_loss_sum(tendency)
        R = LinearRemineralization(parameters.DON_remineralization)(DON)

        frac = one(parameters.DOM_POM_fractionation) - parameters.DOM_POM_fractionation
        return frac * parameters.nitrogen_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""PON tendency assuming fixed stoichiometry (N:C)."""
function PON_default(plankton_syms)
    n_plankton = length(plankton_syms)

    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        PON = tracers.PON(args)

        M = _mortality_loss_sum(tendency, n_plankton)
        g = _grazing_assimilation_loss_sum(tendency)
        R = LinearRemineralization(parameters.PON_remineralization)(PON)

        return parameters.DOM_POM_fractionation * parameters.nitrogen_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""DOP tendency assuming fixed stoichiometry (P:C)."""
function DOP_default(plankton_syms)
    n_plankton = length(plankton_syms)

    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        DOP = tracers.DOP(args)

        M = _mortality_loss_sum(tendency, n_plankton)
        g = _grazing_assimilation_loss_sum(tendency)
        R = LinearRemineralization(parameters.DOP_remineralization)(DOP)

        frac = one(parameters.DOM_POM_fractionation) - parameters.DOM_POM_fractionation
        return frac * parameters.phosphorus_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""POP tendency assuming fixed stoichiometry (P:C)."""
function POP_default(plankton_syms)
    n_plankton = length(plankton_syms)

    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        POP = tracers.POP(args)

        M = _mortality_loss_sum(tendency, n_plankton)
        g = _grazing_assimilation_loss_sum(tendency)
        R = LinearRemineralization(parameters.POP_remineralization)(POP)

        return parameters.DOM_POM_fractionation * parameters.phosphorus_to_carbon * (M + g) - R
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
    requirements = Requirements(;
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
        tendency = TendencyContext(bgc, args)
        parameters = tendency.parameters
        tracers = tendency.tracers

        DIN = tracers.DIN(args)
        PO4 = tracers.PO4(args)
        P = tracers.plankton(args, plankton_idx)
        PAR = tracers.PAR(args)

        growth = TwoNutrientGrowthGeider(
            parameters.maximum_growth_rate[plankton_idx],
            parameters.half_saturation_DIN[plankton_idx],
            parameters.half_saturation_PO4[plankton_idx],
            parameters.photosynthetic_slope[plankton_idx],
            parameters.chlorophyll_to_carbon_ratio[plankton_idx],
        )(
            DIN, PO4, P, PAR
        )

        grazing = _grazing_loss_sum(tendency, P, plankton_idx)

        mort = LinearLoss(parameters.linear_mortality[plankton_idx])(P)

        return growth - grazing - mort
    end

    return CompiledEquation(f, requirements)
end

"""Zooplankton tendency with preferential grazing gain."""
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
