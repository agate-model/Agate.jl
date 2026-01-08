"""Tracer tendency expressions for the NiPiZD biogeochemical model.

This module constructs `Expr` objects that define tracer tendencies for Oceananigans
continuous-form biogeochemistry models.

All generated expressions are allocation-free and suitable for GPU compilation.
"""

module Tracers

export phytoplankton_default,
    phytoplankton_geider_light,
    zooplankton_default,
    nutrient_default,
    nutrient_geider_light,
    detritus_default

"""Return a sum expression for a collection of terms."""
function _sum_expr(terms)
    isempty(terms) && return :(zero(t))

    s = terms[1]
    for i in 2:length(terms)
        s = :($s + $(terms[i]))
    end
    return s
end

"""Build a sum of linear losses over all plankton."""
function _linear_loss_sum(plankton_syms::AbstractVector{Symbol})
    terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(terms, :(linear_loss($sym, linear_mortality[$i])))
    end
    return _sum_expr(terms)
end

"""Build a sum of quadratic losses over all plankton."""
function _quadratic_loss_sum(plankton_syms::AbstractVector{Symbol})
    terms = Expr[]
    for (i, sym) in enumerate(plankton_syms)
        push!(terms, :(quadratic_loss($sym, quadratic_mortality[$i])))
    end
    return _sum_expr(terms)
end

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
    return _sum_expr(terms)
end

"""Build a sum of Geider-style photosynthetic growth over all plankton.

!!! tip
    Nutrients and nutrient half-saturation vectors are passed as single-element tuples e.g. (N,).
""" 
function _photosynthetic_growth_geider_sum(plankton_syms::AbstractVector{Symbol})
    plankton_tuple = Expr(:tuple, plankton_syms...)  # (Z1, Z2, ..., P1, P2, ...)
    return :(photosynthetic_growth_geider_light(
        (N,),                          # nutrients
        (nutrient_half_saturation,),   # per-plankton half-sat vector(s)
        $plankton_tuple,               # per-plankton concentrations
        PAR,
        maximum_growth_rate,
        photosynthetic_slope,
        chlorophyll_to_carbon_ratio,
    ))
end

"""Build the preferential predation loss of a prey size class to all predators."""
function _predation_loss_sum(
    prey_sym::Symbol, prey_idx::Int, plankton_syms::AbstractVector{Symbol}
)
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

"""Build the preferential predation gain of a predator size class from all prey."""
function _predation_gain_sum(
    predator_sym::Symbol, predator_idx::Int, plankton_syms::AbstractVector{Symbol}
)
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

"""Build the total assimilation loss from all predator-prey pairs."""
function _predation_assimilation_loss_sum(plankton_syms::AbstractVector{Symbol})
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

"""
    nutrient_default(plankton_syms)

Nutrient tendency for a single dissolved inorganic nutrient `N`.
"""
function nutrient_default(plankton_syms::AbstractVector{Symbol})
    linear_sum = _linear_loss_sum(plankton_syms)
    quadratic_sum = _quadratic_loss_sum(plankton_syms)
    growth_sum = _photosynthetic_growth_sum(plankton_syms)

    return :(
        mortality_export_fraction * ($linear_sum) +
        mortality_export_fraction * ($quadratic_sum) +
        remineralization_idealized(D, detritus_remineralization) - ($growth_sum)
    )
end

"""Nutrient tendency using Geider-style light limitation."""
function nutrient_geider_light(plankton_syms::AbstractVector{Symbol})
    linear_sum = _linear_loss_sum(plankton_syms)
    quadratic_sum = _quadratic_loss_sum(plankton_syms)
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
    linear_sum = _linear_loss_sum(plankton_syms)
    quadratic_sum = _quadratic_loss_sum(plankton_syms)
    assimilation_loss_sum = _predation_assimilation_loss_sum(plankton_syms)

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

    grazing_loss = _predation_loss_sum(plankton_sym, plankton_idx, plankton_syms)

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

    grazing_loss = _predation_loss_sum(plankton_sym, plankton_idx, plankton_syms)
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
    gain_sum = _predation_gain_sum(plankton_sym, plankton_idx, plankton_syms)

    linear = :(linear_loss($plankton_sym, linear_mortality[$plankton_idx]))
    quadratic = :(quadratic_loss($plankton_sym, quadratic_mortality[$plankton_idx]))

    return :(($gain_sum) - $linear - $quadratic)
end

end # module
