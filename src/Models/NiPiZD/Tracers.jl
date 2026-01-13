"""Tracer tendency equations for the NiPiZD biogeochemical model.

All builders in this module return `Agate.Library.Equations.Equation`.
The symbolic API is used only at construction time; kernels remain plain `Expr`
operating on numeric arrays and scalars.
"""

module Tracers

using Agate.Library.Equations: Equation, Σ, bgc_param

using Agate.Library.Mortality: linear_loss, quadratic_loss, linear_loss_sum, quadratic_loss_sum
using Agate.Library.Predation: grazing_loss, grazing_gain, grazing_assimilation_loss
using Agate.Library.Photosynthesis:
    growth_single_nutrient,
    growth_single_nutrient_comm,
    growth_single_nutrient_geider,
    growth_single_nutrient_geider_comm

export phytoplankton_default,
    phytoplankton_geider_light,
    zooplankton_default,
    nutrient_default,
    nutrient_geider_light,
    detritus_default

# --- internal helpers ---------------------------------------------------------

_growth_sum(plankton_syms) = Σ(plankton_syms) do sym, i
    growth_single_nutrient_comm(:N, sym, :PAR, i)
end

_growth_sum_geider(plankton_syms) = Σ(plankton_syms) do sym, i
    growth_single_nutrient_geider_comm(:N, sym, :PAR, i)
end

_remineralization_term() = bgc_param(:detritus_remineralization) * :D

# --- biogeochemical tracers ---------------------------------------------------

"""Nutrient tendency for a single dissolved inorganic nutrient `N`."""
function nutrient_default(plankton_syms)
    linear_sum = linear_loss_sum(plankton_syms)
    quadratic_sum = quadratic_loss_sum(plankton_syms)
    growth_sum = _growth_sum(plankton_syms)

    export_frac = bgc_param(:mortality_export_fraction)
    remin = _remineralization_term()

    return Equation(export_frac * linear_sum + export_frac * quadratic_sum + remin - growth_sum)
end

"""Nutrient tendency using Geider-style light limitation."""
function nutrient_geider_light(plankton_syms)
    linear_sum = linear_loss_sum(plankton_syms)
    quadratic_sum = quadratic_loss_sum(plankton_syms)
    growth_sum = _growth_sum_geider(plankton_syms)

    export_frac = bgc_param(:mortality_export_fraction)
    remin = _remineralization_term()

    return Equation(export_frac * linear_sum + export_frac * quadratic_sum + remin - growth_sum)
end

"""Detritus tendency for a single detrital pool `D`."""
function detritus_default(plankton_syms)
    linear_sum = linear_loss_sum(plankton_syms)
    quadratic_sum = quadratic_loss_sum(plankton_syms)
    assimilation_loss_sum = grazing_assimilation_loss(plankton_syms)

    export_frac = bgc_param(:mortality_export_fraction)
    remin = _remineralization_term()

    return Equation((1 - export_frac) * linear_sum + assimilation_loss_sum + (1 - export_frac) * quadratic_sum - remin)
end

# --- plankton tracers ---------------------------------------------------------

"""Phytoplankton tendency with single-nutrient Smith-style light limitation."""
function phytoplankton_default(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    growth = growth_single_nutrient(:N, plankton_sym, :PAR, plankton_idx)
    grazing = grazing_loss(plankton_sym, plankton_idx, plankton_syms)
    mort = linear_loss(plankton_sym, plankton_idx)
    return Equation(growth - grazing - mort)
end

"""Phytoplankton tendency using Geider-style light limitation."""
function phytoplankton_geider_light(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    growth = growth_single_nutrient_geider(:N, plankton_sym, :PAR, plankton_idx)
    grazing = grazing_loss(plankton_sym, plankton_idx, plankton_syms)
    mort = linear_loss(plankton_sym, plankton_idx)
    return Equation(growth - grazing - mort)
end

"""Zooplankton tendency with preferential feeding."""
function zooplankton_default(plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    gain = grazing_gain(plankton_sym, plankton_idx, plankton_syms)
    lin = linear_loss(plankton_sym, plankton_idx)
    quad = quadratic_loss(plankton_sym, plankton_idx)
    return Equation(gain - lin - quad)
end

end # module
