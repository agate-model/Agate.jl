"""Tracer tendency expressions for the NiPiZD biogeochemical model.

This module constructs `Expr` objects that define tracer tendencies for Oceananigans
continuous-form biogeochemistry models.

All generated expressions are allocation-free and suitable for GPU compilation.
"""

module Tracers

using Agate.Utils: @register_dynamics, sum_expr
using Agate.Library.Mortality: linear_loss_sum, quadratic_loss_sum
using Agate.Library.Predation: predation_loss_sum, predation_gain_sum, predation_assimilation_loss_sum

export phytoplankton_default,
    phytoplankton_geider_light,
    zooplankton_default,
    nutrient_default,
    nutrient_geider_light,
    detritus_default

"""Build a sum of photosynthetic growth over all plankton."""
function _photosynthetic_growth_sum(plankton_syms::AbstractVector{Symbol})
    terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(
            terms,
            :(photosynthetic_growth_single_nutrient(
                N,
                $sym,
                PAR,
                maximum_growth_rate[$i],
                nutrient_half_saturation[$i],
                alpha[$i],
            )),
        )
    end
    return sum_expr(terms)
end

"""Build a sum of Geider-style photosynthetic growth over all plankton."""
function _photosynthetic_growth_geider_sum(plankton_syms::AbstractVector{Symbol})
    terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(
            terms,
            :(photosynthetic_growth_single_nutrient_geider_light(
                N,
                $sym,
                PAR,
                maximum_growth_rate[$i],
                nutrient_half_saturation[$i],
                photosynthetic_slope[$i],
                chlorophyll_to_carbon_ratio[$i],
            )),
        )
    end
    return sum_expr(terms)
end

"""
    nutrient_default(plankton_syms)

Nutrient tendency for a single dissolved inorganic nutrient `N`.
"""
function nutrient_default(plankton_syms::AbstractVector{Symbol})
    linear_sum = linear_loss_sum(plankton_syms)
    quadratic_sum = quadratic_loss_sum(plankton_syms)
    growth_sum = _photosynthetic_growth_sum(plankton_syms)

    return :(
        mortality_export_fraction * ($linear_sum) +
        mortality_export_fraction * ($quadratic_sum) +
        remineralization_idealized(D, detritus_remineralization) - ($growth_sum)
    )
end

"""Nutrient tendency using Geider-style light limitation."""
function nutrient_geider_light(plankton_syms::AbstractVector{Symbol})
    linear_sum = linear_loss_sum(plankton_syms)
    quadratic_sum = quadratic_loss_sum(plankton_syms)
    growth_sum = _photosynthetic_growth_geider_sum(plankton_syms)

    return :(
        mortality_export_fraction * ($linear_sum) +
        mortality_export_fraction * ($quadratic_sum) +
        remineralization_idealized(D, detritus_remineralization) - ($growth_sum)
    )
end

"""
    detritus_default(plankton_syms)

Detritus tendency for a single detrital pool `D`.
"""
function detritus_default(plankton_syms::AbstractVector{Symbol})
    linear_sum = linear_loss_sum(plankton_syms)
    quadratic_sum = quadratic_loss_sum(plankton_syms)
    assimilation_loss_sum = predation_assimilation_loss_sum(plankton_syms)

    return :(
        (1 - mortality_export_fraction) * ($linear_sum) +
        ($assimilation_loss_sum) +
        (1 - mortality_export_fraction) * ($quadratic_sum) -
        remineralization_idealized(D, detritus_remineralization)
    )
end

"""
    phytoplankton_default(plankton_syms, plankton_sym, plankton_idx)

Phytoplankton tendency with single-nutrient Smith-style light limitation.
"""
function phytoplankton_default(
    plankton_syms::AbstractVector{Symbol}, plankton_sym::Symbol, plankton_idx::Int
)
    growth = :(photosynthetic_growth_single_nutrient(
        N,
        $plankton_sym,
        PAR,
        maximum_growth_rate[$plankton_idx],
        nutrient_half_saturation[$plankton_idx],
        alpha[$plankton_idx],
    ))

    grazing_loss = predation_loss_sum(plankton_syms, plankton_sym, plankton_idx)

    mortality_loss = :(linear_loss($plankton_sym, linear_mortality[$plankton_idx]))

    return :($growth - ($grazing_loss) - $mortality_loss)
end

"""Phytoplankton tendency using Geider-style light limitation."""
function phytoplankton_geider_light(
    plankton_syms::AbstractVector{Symbol}, plankton_sym::Symbol, plankton_idx::Int
)
    growth = :(photosynthetic_growth_single_nutrient_geider_light(
        N,
        $plankton_sym,
        PAR,
        maximum_growth_rate[$plankton_idx],
        nutrient_half_saturation[$plankton_idx],
        photosynthetic_slope[$plankton_idx],
        chlorophyll_to_carbon_ratio[$plankton_idx],
    ))

    grazing_loss = predation_loss_sum(plankton_syms, plankton_sym, plankton_idx)
    mortality_loss = :(linear_loss($plankton_sym, linear_mortality[$plankton_idx]))

    return :($growth - ($grazing_loss) - $mortality_loss)
end

"""
    zooplankton_default(plankton_syms, plankton_sym, plankton_idx)

Zooplankton tendency with preferential feeding.
"""
function zooplankton_default(
    plankton_syms::AbstractVector{Symbol}, plankton_sym::Symbol, plankton_idx::Int
)
    gain_sum = predation_gain_sum(plankton_syms, plankton_sym, plankton_idx)

    linear = :(linear_loss($plankton_sym, linear_mortality[$plankton_idx]))
    quadratic = :(quadratic_loss($plankton_sym, quadratic_mortality[$plankton_idx]))

    return :(($gain_sum) - $linear - $quadratic)
end

# -----------------------------------------------------------------------------
# Dynamics parameter registry
# -----------------------------------------------------------------------------

@register_dynamics nutrient_default (
    :linear_mortality,
    :quadratic_mortality,
    :maximum_growth_rate,
    :nutrient_half_saturation,
    :alpha,
    :detritus_remineralization,
    :mortality_export_fraction,
)

@register_dynamics nutrient_geider_light (
    :linear_mortality,
    :quadratic_mortality,
    :maximum_growth_rate,
    :nutrient_half_saturation,
    :photosynthetic_slope,
    :chlorophyll_to_carbon_ratio,
    :detritus_remineralization,
    :mortality_export_fraction,
)

@register_dynamics detritus_default (
    :linear_mortality,
    :quadratic_mortality,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :assimilation_efficiency_matrix,
    :detritus_remineralization,
    :mortality_export_fraction,
)

@register_dynamics phytoplankton_default (
    :maximum_growth_rate,
    :nutrient_half_saturation,
    :alpha,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :linear_mortality,
)

@register_dynamics phytoplankton_geider_light (
    :maximum_growth_rate,
    :nutrient_half_saturation,
    :photosynthetic_slope,
    :chlorophyll_to_carbon_ratio,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :linear_mortality,
)

@register_dynamics zooplankton_default (
    :assimilation_efficiency_matrix,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :linear_mortality,
    :quadratic_mortality,
)

end # module
