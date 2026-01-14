"""Tracer tendency equations for the simplified DARWIN-like elemental cycling model.

All builders in this module return `Agate.Library.Equations.Equation`.
The symbolic API is used only at construction time; kernels remain plain `Expr`
operating on numeric arrays and scalars.
"""

module Tracers

using Agate.Library.Equations: Equation, Σ

using Agate.Library.Mortality: linear_loss, quadratic_loss, linear_loss_sum, quadratic_loss_sum
using Agate.Library.Predation: grazing_loss, grazing_gain, grazing_assimilation_loss
using Agate.Library.Photosynthesis: growth_two_nutrients_geider
using Agate.Library.Remineralization: remineralization_flux

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
_growth_sum(plankton_syms) = Σ(plankton_syms) do sym, i
    growth_two_nutrients_geider(:DIN, :PO4, sym, :PAR, i)
end

"""Base loss term contributing to organic matter: linear + quadratic + sloppy-feeding losses."""
_base_loss(plankton_syms) = linear_loss_sum(plankton_syms) +
                            grazing_assimilation_loss(plankton_syms) +
                            quadratic_loss_sum(plankton_syms)

# ----------------------------------------------------------------------------
# Inorganic tracers
# ----------------------------------------------------------------------------

"""DIC tendency with Geider-style growth (carbon units)."""
function DIC_geider_light(plankton_syms)
    growth = _growth_sum(plankton_syms)
    return Equation(
        remineralization_flux(:DOC, :DOC_remineralization) +
        remineralization_flux(:POC, :POC_remineralization) -
        growth,
    )
end

"""DIN tendency assuming fixed stoichiometry (N:C)."""
function DIN_geider_light(plankton_syms)
    growth = _growth_sum(plankton_syms)
    return Equation(
        remineralization_flux(:DON, :DON_remineralization) +
        remineralization_flux(:PON, :PON_remineralization) -
        nitrogen_to_carbon * growth,
    )
end

"""PO4 tendency assuming fixed stoichiometry (P:C)."""
function PO4_geider_light(plankton_syms)
    growth = _growth_sum(plankton_syms)
    return Equation(
        remineralization_flux(:DOP, :DOP_remineralization) +
        remineralization_flux(:POP, :POP_remineralization) -
        phosphorus_to_carbon * growth,
    )
end

# ----------------------------------------------------------------------------
# Organic matter tracers
# ----------------------------------------------------------------------------

"""DOC tendency from plankton losses and remineralization."""
function DOC_default(plankton_syms)
    base = _base_loss(plankton_syms)
    return Equation(
        (1 - DOM_POM_fractionation) * base -
        remineralization_flux(:DOC, :DOC_remineralization),
    )
end

"""POC tendency from plankton losses and remineralization."""
function POC_default(plankton_syms)
    base = _base_loss(plankton_syms)
    return Equation(
        DOM_POM_fractionation * base -
        remineralization_flux(:POC, :POC_remineralization),
    )
end

"""DON tendency assuming fixed stoichiometry (N:C)."""
function DON_default(plankton_syms)
    base = _base_loss(plankton_syms)
    return Equation(
        (1 - DOM_POM_fractionation) * nitrogen_to_carbon * base -
        remineralization_flux(:DON, :DON_remineralization),
    )
end

"""PON tendency assuming fixed stoichiometry (N:C)."""
function PON_default(plankton_syms)
    base = _base_loss(plankton_syms)
    return Equation(
        DOM_POM_fractionation * nitrogen_to_carbon * base -
        remineralization_flux(:PON, :PON_remineralization),
    )
end

"""DOP tendency assuming fixed stoichiometry (P:C)."""
function DOP_default(plankton_syms)
    base = _base_loss(plankton_syms)
    return Equation(
        (1 - DOM_POM_fractionation) * phosphorus_to_carbon * base -
        remineralization_flux(:DOP, :DOP_remineralization),
    )
end

"""POP tendency assuming fixed stoichiometry (P:C)."""
function POP_default(plankton_syms)
    base = _base_loss(plankton_syms)
    return Equation(
        DOM_POM_fractionation * phosphorus_to_carbon * base -
        remineralization_flux(:POP, :POP_remineralization),
    )
end

# ----------------------------------------------------------------------------
# Plankton tracers
# ----------------------------------------------------------------------------

"""Phytoplankton tendency with Geider-style, two-nutrient growth."""
function phytoplankton_growth_two_nutrients_geider_light(
    plankton_syms,
    plankton_sym::Symbol,
    plankton_idx::Int,
)
    growth = growth_two_nutrients_geider(:DIN, :PO4, plankton_sym, :PAR, plankton_idx)
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
