module Tracers

using Agate.Utils: @register_dynamics

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

# -----------------------------------------------------------------------------
# Internal helpers to build allocation-free expressions
# -----------------------------------------------------------------------------

@inline function _sum_expr(terms::Vector{Expr})
    if isempty(terms)
        return :(0)
    end
    out = terms[1]
    @inbounds for i in 2:length(terms)
        out = :($out + $(terms[i]))
    end
    return out
end

@inline function _linear_loss_sum(plankton_syms)
    terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(terms, :(linear_loss($sym, linear_mortality[$i])))
    end
    return _sum_expr(terms)
end

@inline function _quadratic_loss_sum(plankton_syms)
    terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(terms, :(quadratic_loss($sym, quadratic_mortality[$i])))
    end
    return _sum_expr(terms)
end

@inline function _photosynthetic_growth_sum(plankton_syms)
    terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(
            terms,
            :(photosynthetic_growth_two_nutrients_geider_light(
                DIN,
                PO4,
                $sym,
                PAR,
                maximum_growth_rate[$i],
                half_saturation_DIN[$i],
                half_saturation_PO4[$i],
                photosynthetic_slope[$i],
                chlorophyll_to_carbon_ratio[$i],
            )),
        )
    end
    return _sum_expr(terms)
end

@inline function _predation_loss_sum(plankton_syms, prey_sym, prey_idx)
    terms = Expr[]
    for (pred_idx, pred_sym) in enumerate(plankton_syms)
        push!(
            terms,
            :(predation_loss_preferential(
                $prey_sym,
                $pred_sym,
                maximum_predation_rate[$pred_idx],
                holling_half_saturation[$pred_idx],
                palatability_matrix[$pred_idx, $prey_idx],
            )),
        )
    end
    return _sum_expr(terms)
end

@inline function _predation_gain_sum(plankton_syms, predator_sym, predator_idx)
    terms = Expr[]
    for (prey_idx, prey_sym) in enumerate(plankton_syms)
        push!(
            terms,
            :(predation_gain_preferential(
                $prey_sym,
                $predator_sym,
                assimilation_efficiency_matrix[$predator_idx, $prey_idx],
                maximum_predation_rate[$predator_idx],
                holling_half_saturation[$predator_idx],
                palatability_matrix[$predator_idx, $prey_idx],
            )),
        )
    end
    return _sum_expr(terms)
end

@inline function _predation_assimilation_loss_sum(plankton_syms)
    terms = Expr[]
    for (pred_idx, pred_sym) in enumerate(plankton_syms)
        for (prey_idx, prey_sym) in enumerate(plankton_syms)
            push!(
                terms,
                :(predation_assimilation_loss_preferential(
                    $prey_sym,
                    $pred_sym,
                    assimilation_efficiency_matrix[$pred_idx, $prey_idx],
                    maximum_predation_rate[$pred_idx],
                    holling_half_saturation[$pred_idx],
                    palatability_matrix[$pred_idx, $prey_idx],
                )),
            )
        end
    end
    return _sum_expr(terms)
end

# -----------------------------------------------------------------------------
# Inorganic tracers
# -----------------------------------------------------------------------------

"""
    DIC_geider_light(plankton_syms)

Build expression representing the evolution of DIC over time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function DIC_geider_light(plankton_syms)
    growth = _photosynthetic_growth_sum(plankton_syms)
    return :(
        remineralization_idealized(DOC, DOC_remineralization) +
        remineralization_idealized(POC, POC_remineralization) - ($growth)
    )
end

"""
    DIN_geider_light(plankton_syms)

Build expression representing the evolution of DIN over time assuming fixed stoichiometry.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function DIN_geider_light(plankton_syms)
    growth = _photosynthetic_growth_sum(plankton_syms)
    return :(
        remineralization_idealized(DON, DON_remineralization) +
        remineralization_idealized(PON, PON_remineralization) -
        nitrogen_to_carbon * ($growth)
    )
end

"""
    PO4_geider_light(plankton_syms)

Build expression representing the evolution of PO4 over time assuming fixed stoichiometry.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function PO4_geider_light(plankton_syms)
    growth = _photosynthetic_growth_sum(plankton_syms)
    return :(
        remineralization_idealized(DOP, DOP_remineralization) +
        remineralization_idealized(POP, POP_remineralization) -
        phosphorus_to_carbon * ($growth)
    )
end

# -----------------------------------------------------------------------------
# Organic matter tracers
# -----------------------------------------------------------------------------

"""
    DOC_default(plankton_syms)

Build expression for a simplified DOC function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function DOC_default(plankton_syms)
    base = :(
        $(_linear_loss_sum(plankton_syms)) +
        $(_predation_assimilation_loss_sum(plankton_syms)) +
        $(_quadratic_loss_sum(plankton_syms))
    )
    return :(
        (1 - DOM_POM_fractionation) * ($base) -
        remineralization_idealized(DOC, DOC_remineralization)
    )
end

"""
    POC_default(plankton_syms)

Build expression for a simplified POC function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function POC_default(plankton_syms)
    base = :(
        $(_linear_loss_sum(plankton_syms)) +
        $(_predation_assimilation_loss_sum(plankton_syms)) +
        $(_quadratic_loss_sum(plankton_syms))
    )
    return :(
        DOM_POM_fractionation * ($base) -
        remineralization_idealized(POC, POC_remineralization)
    )
end

"""
    DON_default(plankton_syms)

Build expression for a simplified DON function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function DON_default(plankton_syms)
    base = :(
        $(_linear_loss_sum(plankton_syms)) +
        $(_predation_assimilation_loss_sum(plankton_syms)) +
        $(_quadratic_loss_sum(plankton_syms))
    )
    return :(
        (1 - DOM_POM_fractionation) * nitrogen_to_carbon * ($base) -
        remineralization_idealized(DON, DON_remineralization)
    )
end

"""
    PON_default(plankton_syms)

Build expression for a simplified PON function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function PON_default(plankton_syms)
    base = :(
        $(_linear_loss_sum(plankton_syms)) +
        $(_predation_assimilation_loss_sum(plankton_syms)) +
        $(_quadratic_loss_sum(plankton_syms))
    )
    return :(
        DOM_POM_fractionation * nitrogen_to_carbon * ($base) -
        remineralization_idealized(PON, PON_remineralization)
    )
end

"""
    DOP_default(plankton_syms)

Build expression for a simplified DOP function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function DOP_default(plankton_syms)
    base = :(
        $(_linear_loss_sum(plankton_syms)) +
        $(_predation_assimilation_loss_sum(plankton_syms)) +
        $(_quadratic_loss_sum(plankton_syms))
    )
    return :(
        (1 - DOM_POM_fractionation) * phosphorus_to_carbon * ($base) -
        remineralization_idealized(DOP, DOP_remineralization)
    )
end

"""
    POP_default(plankton_syms)

Build expression for a simplified POP function of time.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
"""
function POP_default(plankton_syms)
    base = :(
        $(_linear_loss_sum(plankton_syms)) +
        $(_predation_assimilation_loss_sum(plankton_syms)) +
        $(_quadratic_loss_sum(plankton_syms))
    )
    return :(
        DOM_POM_fractionation * phosphorus_to_carbon * ($base) -
        remineralization_idealized(POP, POP_remineralization)
    )
end

# -----------------------------------------------------------------------------
# Plankton tracers
# -----------------------------------------------------------------------------

"""
    phytoplankton_growth_two_nutrients_geider_light(plankton_syms, plankton_sym, plankton_idx)

Build expression for a simplified phytoplankton growth function.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
- `plankton_sym`: name of the phytoplankton for which we are returning the expression, passed
    as a Symbol (e.g., :P1).
- `plankton_idx`: the index at which this plankton's values are stored in all parameter Arrays
"""
function phytoplankton_growth_two_nutrients_geider_light(
    plankton_syms, plankton_sym, plankton_idx
)
    growth = :(photosynthetic_growth_two_nutrients_geider_light(
        DIN,
        PO4,
        $plankton_sym,
        PAR,
        maximum_growth_rate[$plankton_idx],
        half_saturation_DIN[$plankton_idx],
        half_saturation_PO4[$plankton_idx],
        photosynthetic_slope[$plankton_idx],
        chlorophyll_to_carbon_ratio[$plankton_idx],
    ))

    grazing_loss = _predation_loss_sum(plankton_syms, plankton_sym, plankton_idx)
    lin = :(linear_loss($plankton_sym, linear_mortality[$plankton_idx]))

    return :($growth - ($grazing_loss) - $lin)
end

"""
    zooplankton_default(plankton_syms, plankton_sym, plankton_idx)

Build expression for simplified zooplankton growth function.

The functions used in the expression are all within the Agate.Library, see their docstring
for overview. All arguments in the functions are either an Array or a Float. The Arrays have
to be of same length for correct indexing to work (and arranged in the same plankton order).

# Arguments
- `plankton_syms`: names of all the plankton in the ecosystem expressed as Symbols, e.g.:
    `[:Z1, :Z2, :P1, :P2]`, arranged in the same order as all the parameter Arrays
- `plankton_sym`: name of the zooplankton for which we are returning the expression, passed
    as a Symbol (e.g., :Z1).
- `plankton_idx`: the index at which this plankton's values are stored in all parameter Arrays
"""
function zooplankton_default(plankton_syms, plankton_sym, plankton_idx)
    gain = _predation_gain_sum(plankton_syms, plankton_sym, plankton_idx)
    lin = :(linear_loss($plankton_sym, linear_mortality[$plankton_idx]))
    quad = :(quadratic_loss($plankton_sym, quadratic_mortality[$plankton_idx]))
    return :(($gain) - $lin - $quad)
end

# -----------------------------------------------------------------------------
# Parameter registry
# -----------------------------------------------------------------------------

@register_dynamics DIC_geider_light (
    :DOC_remineralization,
    :POC_remineralization,
    :maximum_growth_rate,
    :half_saturation_DIN,
    :half_saturation_PO4,
    :photosynthetic_slope,
    :chlorophyll_to_carbon_ratio,
)

@register_dynamics DIN_geider_light (
    :DON_remineralization,
    :PON_remineralization,
    :nitrogen_to_carbon,
    :maximum_growth_rate,
    :half_saturation_DIN,
    :half_saturation_PO4,
    :photosynthetic_slope,
    :chlorophyll_to_carbon_ratio,
)

@register_dynamics PO4_geider_light (
    :DOP_remineralization,
    :POP_remineralization,
    :phosphorus_to_carbon,
    :maximum_growth_rate,
    :half_saturation_DIN,
    :half_saturation_PO4,
    :photosynthetic_slope,
    :chlorophyll_to_carbon_ratio,
)

@register_dynamics DOC_default (
    :DOM_POM_fractionation,
    :DOC_remineralization,
    :linear_mortality,
    :quadratic_mortality,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :assimilation_efficiency_matrix,
)

@register_dynamics POC_default (
    :DOM_POM_fractionation,
    :POC_remineralization,
    :linear_mortality,
    :quadratic_mortality,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :assimilation_efficiency_matrix,
)

@register_dynamics DON_default (
    :DOM_POM_fractionation,
    :DON_remineralization,
    :nitrogen_to_carbon,
    :linear_mortality,
    :quadratic_mortality,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :assimilation_efficiency_matrix,
)

@register_dynamics PON_default (
    :DOM_POM_fractionation,
    :PON_remineralization,
    :nitrogen_to_carbon,
    :linear_mortality,
    :quadratic_mortality,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :assimilation_efficiency_matrix,
)

@register_dynamics DOP_default (
    :DOM_POM_fractionation,
    :DOP_remineralization,
    :phosphorus_to_carbon,
    :linear_mortality,
    :quadratic_mortality,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :assimilation_efficiency_matrix,
)

@register_dynamics POP_default (
    :DOM_POM_fractionation,
    :POP_remineralization,
    :phosphorus_to_carbon,
    :linear_mortality,
    :quadratic_mortality,
    :maximum_predation_rate,
    :holling_half_saturation,
    :palatability_matrix,
    :assimilation_efficiency_matrix,
)

@register_dynamics phytoplankton_growth_two_nutrients_geider_light (
    :maximum_growth_rate,
    :half_saturation_DIN,
    :half_saturation_PO4,
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
