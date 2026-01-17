"""Tracer tendency equations for the simplified DARWIN-like elemental cycling model.

All builders in this module return `Agate.Equations.Equation`.
The symbolic API is used only at construction time; kernels remain plain `Expr`
operating on numeric arrays and scalars.

All equation builders accept a first argument `PV`, a namespace of construction-time
parameter placeholders provided by `construct`.
"""

module Tracers

using ....Equations: Equation, sum_over

using ....Library.Mortality: linear_loss, quadratic_loss, linear_loss_sum, quadratic_loss_sum
using ....Library.Predation: grazing_loss, grazing_gain, grazing_assimilation_loss
using ....Library.Photosynthesis: growth_two_nutrients_geider
using ....Library.Remineralization: remineralization_flux

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

# ----------------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------------

"""Sum photosynthetic growth across all plankton."""
_growth_sum(PV, plankton_syms) = sum_over(plankton_syms) do sym, i
    growth_two_nutrients_geider(PV, :DIN, :PO4, sym, :PAR, i)
end

"""Base loss term contributing to organic matter: linear + quadratic + sloppy-feeding losses."""
_base_loss(PV, plankton_syms) = linear_loss_sum(PV, plankton_syms) +
                               grazing_assimilation_loss(PV, plankton_syms) +
                               quadratic_loss_sum(PV, plankton_syms)

# ----------------------------------------------------------------------------
# Inorganic tracers
# ----------------------------------------------------------------------------

"""DIC tendency with Geider-style growth (carbon units)."""
function DIC_geider_light(PV, plankton_syms)
    growth = _growth_sum(PV, plankton_syms)
    return Equation(
        remineralization_flux(PV, :DOC, :DOC_remineralization) +
        remineralization_flux(PV, :POC, :POC_remineralization) -
        growth,
    )
end

"""DIN tendency assuming fixed stoichiometry (N:C)."""
function DIN_geider_light(PV, plankton_syms)
    growth = _growth_sum(PV, plankton_syms)
    return Equation(
        remineralization_flux(PV, :DON, :DON_remineralization) +
        remineralization_flux(PV, :PON, :PON_remineralization) -
        PV.nitrogen_to_carbon * growth,
    )
end

"""PO4 tendency assuming fixed stoichiometry (P:C)."""
function PO4_geider_light(PV, plankton_syms)
    growth = _growth_sum(PV, plankton_syms)
    return Equation(
        remineralization_flux(PV, :DOP, :DOP_remineralization) +
        remineralization_flux(PV, :POP, :POP_remineralization) -
        PV.phosphorus_to_carbon * growth,
    )
end

# ----------------------------------------------------------------------------
# Organic matter tracers
# ----------------------------------------------------------------------------

"""DOC tendency from plankton losses and remineralization."""
function DOC_default(PV, plankton_syms)
    base = _base_loss(PV, plankton_syms)
    return Equation(
        (1 - PV.DOM_POM_fractionation) * base -
        remineralization_flux(PV, :DOC, :DOC_remineralization),
    )
end

"""POC tendency from plankton losses and remineralization."""
function POC_default(PV, plankton_syms)
    base = _base_loss(PV, plankton_syms)
    return Equation(
        PV.DOM_POM_fractionation * base -
        remineralization_flux(PV, :POC, :POC_remineralization),
    )
end

"""DON tendency assuming fixed stoichiometry (N:C)."""
function DON_default(PV, plankton_syms)
    base = _base_loss(PV, plankton_syms)
    return Equation(
        (1 - PV.DOM_POM_fractionation) * PV.nitrogen_to_carbon * base -
        remineralization_flux(PV, :DON, :DON_remineralization),
    )
end

"""PON tendency assuming fixed stoichiometry (N:C)."""
function PON_default(PV, plankton_syms)
    base = _base_loss(PV, plankton_syms)
    return Equation(
        PV.DOM_POM_fractionation * PV.nitrogen_to_carbon * base -
        remineralization_flux(PV, :PON, :PON_remineralization),
    )
end

"""DOP tendency assuming fixed stoichiometry (P:C)."""
function DOP_default(PV, plankton_syms)
    base = _base_loss(PV, plankton_syms)
    return Equation(
        (1 - PV.DOM_POM_fractionation) * PV.phosphorus_to_carbon * base -
        remineralization_flux(PV, :DOP, :DOP_remineralization),
    )
end

"""POP tendency assuming fixed stoichiometry (P:C)."""
function POP_default(PV, plankton_syms)
    base = _base_loss(PV, plankton_syms)
    return Equation(
        PV.DOM_POM_fractionation * PV.phosphorus_to_carbon * base -
        remineralization_flux(PV, :POP, :POP_remineralization),
    )
end

# ----------------------------------------------------------------------------
# Plankton tracers
# ----------------------------------------------------------------------------

"""Phytoplankton tendency with Geider-style, two-nutrient growth."""
function phytoplankton_growth_two_nutrients_geider_light(
    PV,
    plankton_syms,
    plankton_sym::Symbol,
    plankton_idx::Int,
)
    growth = growth_two_nutrients_geider(PV, :DIN, :PO4, plankton_sym, :PAR, plankton_idx)
    grazing = grazing_loss(PV, plankton_sym, plankton_idx, plankton_syms)
    mort = linear_loss(PV, plankton_sym, plankton_idx)
    return Equation(growth - grazing - mort)
end

"""Zooplankton tendency with preferential feeding."""
function zooplankton_default(PV, plankton_syms, plankton_sym::Symbol, plankton_idx::Int)
    gain = grazing_gain(PV, plankton_sym, plankton_idx, plankton_syms)
    lin = linear_loss(PV, plankton_sym, plankton_idx)
    quad = quadratic_loss(PV, plankton_sym, plankton_idx)
    return Equation(gain - lin - quad)
end

end # module
