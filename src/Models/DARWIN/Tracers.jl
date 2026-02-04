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

using ....Library.Mortality: linear_loss, quadratic_loss
using ....Library.Photosynthesis: geider_two_nutrient_growth
using ....Library.Remineralization: linear_remineralization

using ....Utils: tendency_views

using ...Sums:
    grazing_unassimilated_loss_sum,
    grazing_loss_sum,
    grazing_gain_sum,
    geider_two_nutrient_uptake_sum,
    mortality_loss_sum

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

"""DIC tendency with Geider-style growth (carbon units)."""
function DIC_geider_light(plankton_syms)
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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton

        DIN = vals.DIN
        PO4 = vals.PO4
        DOC = vals.DOC
        POC = vals.POC
        PAR = vals.PAR

        uptake = geider_two_nutrient_uptake_sum(
            plankton,
            DIN,
            PO4,
            PAR,
            parameters.maximum_growth_rate,
            parameters.half_saturation_DIN,
            parameters.half_saturation_PO4,
            parameters.photosynthetic_slope,
            parameters.chlorophyll_to_carbon_ratio,
        )

        dic_remin =
            linear_remineralization(DOC, parameters.DOC_remineralization) +
            linear_remineralization(POC, parameters.POC_remineralization)

        return dic_remin - uptake
    end

    return CompiledEquation(f, requirements)
end

"""DIN tendency assuming fixed stoichiometry (N:C)."""
function DIN_geider_light(plankton_syms)
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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton

        DIN = vals.DIN
        PO4 = vals.PO4
        DON = vals.DON
        PON = vals.PON
        PAR = vals.PAR

        uptake = geider_two_nutrient_uptake_sum(
            plankton,
            DIN,
            PO4,
            PAR,
            parameters.maximum_growth_rate,
            parameters.half_saturation_DIN,
            parameters.half_saturation_PO4,
            parameters.photosynthetic_slope,
            parameters.chlorophyll_to_carbon_ratio,
        )

        din_remin =
            linear_remineralization(DON, parameters.DON_remineralization) +
            linear_remineralization(PON, parameters.PON_remineralization)

        return din_remin - parameters.nitrogen_to_carbon * uptake
    end

    return CompiledEquation(f, requirements)
end

"""PO4 tendency assuming fixed stoichiometry (P:C)."""
function PO4_geider_light(plankton_syms)
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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton

        DIN = vals.DIN
        PO4 = vals.PO4
        DOP = vals.DOP
        POP = vals.POP
        PAR = vals.PAR

        uptake = geider_two_nutrient_uptake_sum(
            plankton,
            DIN,
            PO4,
            PAR,
            parameters.maximum_growth_rate,
            parameters.half_saturation_DIN,
            parameters.half_saturation_PO4,
            parameters.photosynthetic_slope,
            parameters.chlorophyll_to_carbon_ratio,
        )

        po4_remin =
            linear_remineralization(DOP, parameters.DOP_remineralization) +
            linear_remineralization(POP, parameters.POP_remineralization)

        return po4_remin - parameters.phosphorus_to_carbon * uptake
    end

    return CompiledEquation(f, requirements)
end

# --- Organic matter ---------------------------------------------------------

"""DOC tendency from plankton losses and remineralization."""
function DOC_default(plankton_syms)
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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton
        DOC = vals.DOC

        M = mortality_loss_sum(
            plankton,
            parameters.linear_mortality,
            parameters.quadratic_mortality,
            zero(DOC),
        )
        g = grazing_unassimilated_loss_sum(parameters, plankton)
        R = linear_remineralization(DOC, parameters.DOC_remineralization)

        frac = one(parameters.DOM_POM_fractionation) - parameters.DOM_POM_fractionation
        return frac * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""POC tendency from plankton losses and remineralization."""
function POC_default(plankton_syms)
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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton
        POC = vals.POC

        M = mortality_loss_sum(
            plankton,
            parameters.linear_mortality,
            parameters.quadratic_mortality,
            zero(POC),
        )
        g = grazing_unassimilated_loss_sum(parameters, plankton)
        R = linear_remineralization(POC, parameters.POC_remineralization)

        return parameters.DOM_POM_fractionation * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""DON tendency assuming fixed stoichiometry (N:C)."""
function DON_default(plankton_syms)
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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton
        DON = vals.DON

        M = mortality_loss_sum(
            plankton,
            parameters.linear_mortality,
            parameters.quadratic_mortality,
            zero(DON),
        )
        g = grazing_unassimilated_loss_sum(parameters, plankton)
        R = linear_remineralization(DON, parameters.DON_remineralization)

        frac = one(parameters.DOM_POM_fractionation) - parameters.DOM_POM_fractionation
        return frac * parameters.nitrogen_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""PON tendency assuming fixed stoichiometry (N:C)."""
function PON_default(plankton_syms)
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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton
        PON = vals.PON

        M = mortality_loss_sum(
            plankton,
            parameters.linear_mortality,
            parameters.quadratic_mortality,
            zero(PON),
        )
        g = grazing_unassimilated_loss_sum(parameters, plankton)
        R = linear_remineralization(PON, parameters.PON_remineralization)

        return parameters.DOM_POM_fractionation * parameters.nitrogen_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""DOP tendency assuming fixed stoichiometry (P:C)."""
function DOP_default(plankton_syms)
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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton
        DOP = vals.DOP

        M = mortality_loss_sum(
            plankton,
            parameters.linear_mortality,
            parameters.quadratic_mortality,
            zero(DOP),
        )
        g = grazing_unassimilated_loss_sum(parameters, plankton)
        R = linear_remineralization(DOP, parameters.DOP_remineralization)

        frac = one(parameters.DOM_POM_fractionation) - parameters.DOM_POM_fractionation
        return frac * parameters.phosphorus_to_carbon * (M + g) - R
    end

    return CompiledEquation(f, requirements)
end

"""POP tendency assuming fixed stoichiometry (P:C)."""
function POP_default(plankton_syms)
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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton
        POP = vals.POP

        M = mortality_loss_sum(
            plankton,
            parameters.linear_mortality,
            parameters.quadratic_mortality,
            zero(POP),
        )
        g = grazing_unassimilated_loss_sum(parameters, plankton)
        R = linear_remineralization(POP, parameters.POP_remineralization)

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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton

        DIN = vals.DIN
        PO4 = vals.PO4
        PAR = vals.PAR
        P = plankton(plankton_idx)

        μmax = parameters.maximum_growth_rate[plankton_idx]
        KDIN = parameters.half_saturation_DIN[plankton_idx]
        KPO4 = parameters.half_saturation_PO4[plankton_idx]
        α = parameters.photosynthetic_slope[plankton_idx]
        θc = parameters.chlorophyll_to_carbon_ratio[plankton_idx]

        growth = geider_two_nutrient_growth(DIN, PO4, P, PAR, μmax, KDIN, KPO4, α, θc)

        grazing = grazing_loss_sum(parameters, plankton, P, plankton_idx)

        m = parameters.linear_mortality[plankton_idx]
        mort = linear_loss(P, m)

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
        _, parameters, vals = tendency_views(bgc, args)

        plankton = vals.plankton
        Z = plankton(plankton_idx)

        gain = grazing_gain_sum(parameters, plankton, Z, plankton_idx)

        m_lin = parameters.linear_mortality[plankton_idx]
        m_quad = parameters.quadratic_mortality[plankton_idx]

        lin = linear_loss(Z, m_lin)
        quad = quadratic_loss(Z, m_quad)

        return gain - lin - quad
    end

    return CompiledEquation(f, requirements)
end

end # module
