module Predation

using ..ExprUtils: sum_expr

export AssimilationPreyParameters
export AssimilationPredatorParameters
export holling_type_2
export predation_loss_idealized
export predation_gain_idealized
export predation_assimilation_loss_idealized
export predation_loss_preferential
export predation_gain_preferential
export predation_assimilation_loss_preferential
export assimilation_efficiency_emergent_binary
export predation_loss_sum
export predation_gain_sum
export predation_assimilation_loss_sum

"""
    holling_type_2(prey_concentration, prey_half_saturation)

Holling's "type II" functional response as describe in Holling 1959.

!!! formulation
    R / (kᵣ + R)

    where:
    - ``R`` = prey concentration 
    - ``kᵣ`` = prey half saturation constant    

The function is similar to the Monod equation and Michaelis-Menten equation of for enzyme kinetics.
The formulation is characterized by decelerating predation as prey concentrations increase.

# Arguments
- `prey_concentration`: prey density
- `prey_half_saturation`: prey density at which predation is half it's maximum rate
"""
@inline function holling_type_2(prey_concentration, prey_half_saturation)
    return prey_concentration / (prey_half_saturation + prey_concentration)
end

"""
    predation_loss_idealized(P, Z, maximum_grazing_rate, prey_half_saturation)

Estimates the loss rate of P (prey), to Z (predator).

!!! formulation
    gₘₐₓ * γᴾᴿᴱᴰ * Z

    where:
    - gₘₐₓ = maximum grazing rate of the predator
    - γᴾᴿᴱᴰ = `holling_type_2(P², kₚ²)`
    - P = prey concentration
    - kₚ = prey half saturation
    - Z = predator concentration

In this formulation predator-prey interactions are modulated both by their density (Holling type 2)
and the prey-predator palatability.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `maximum_grazing_rate`: maximum grazing rate of the predator
- `prey_half_saturation`: prey density at which predation is half it's maximum rate
"""
@inline function predation_loss_idealized(P, Z, maximum_grazing_rate, prey_half_saturation)
    return maximum_grazing_rate * holling_type_2(P^2, prey_half_saturation^2) * Z
end

"""
    predation_gain_idealized(P, Z, assimilation_efficiency, maximum_grazing_rate, 
        prey_half_saturation)

Estimates the gain rate of Z (predator) feeding on P (prey).
In this formulation predation rate is multiplied by an assimilation efficiency (β) which
represents 'sloppy feeding'.

!!! formulation
    β * ``g``
    
    where:
    - β = assimilation efficiency of prey to the predator
    - ``g`` = `predation_loss_idealized(P, Z, gₘₐₓ, kₚ)`
    - P = phytoplankton concentration
    - Z = zooplankton concentration
    - gₘₐₓ = maximum grazing rate of the predator
    - kₚ = prey half saturation

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `assimilation_efficiency`: assimilation efficiency of the predator
- `maximum_grazing_rate`: maximum grazing rate of the predator
- `prey_half_saturation`: prey half saturation
"""
@inline function predation_gain_idealized(
    P, Z, assimilation_efficiency, maximum_grazing_rate, k_p
)
    return assimilation_efficiency *
           predation_loss_idealized(P, Z, maximum_grazing_rate, k_p)
end

"""
    predation_assimilation_loss_idealized(P, Z, assimilation_efficiency, 
        maximum_grazing_rate, prey_half_saturation)

Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

!!! formulation

    (1 - β) * ``g``
    
    where:
    - β = assimilation efficiency of prey to the predator
    - ``g`` = `predation_loss_idealized(P, Z, gₘₐₓ, kₚ)`
    - P = phytoplankton concentration
    - Z = zooplankton concentration
    - gₘₐₓ = maximum grazing rate of the predator
    - kₚ = prey half saturation

!!! note
    This function differs from `predation_gain_idealized` as it represents the transfer 
    of biomass from the prey to the environment rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `assimilation_efficiency`: assimilation efficiency of prey to the predator
- `maximum_grazing_rate`: maximum grazing rate of the predator
- `prey_half_saturation`: prey half saturation
"""
@inline function predation_assimilation_loss_idealized(
    P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation
)
    return (1 - assimilation_efficiency) *
           predation_loss_idealized(P, Z, maximum_grazing_rate, prey_half_saturation)
end

"""
    predation_loss_preferential(P, Z, maximum_grazing_rate, prey_half_saturation, 
        palatability)

Estimates the loss rate of P (prey), to Z (predator).
In this formulation predator-prey interactions are modulated both by their density (Holling type 2)
and the prey-predator palatability.

!!! formulation
    gₘₐₓ * η * γᴾᴿᴱᴰ * Z

    where:
    - gₘₐₓ = maximum grazing rate
    - η = palatability
    - γᴾᴿᴱᴰ = `holling_type_2(P², kₚ²)`
    - P = phytoplankton concentration
    - kₚ = prey half saturation
    - Z: zooplankton concentration

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `maximum_grazing_rate`: maximum grazing rate
- `prey_half_saturation`: prey density at which predation is half it's maximum rate
- `palatability`: the likelihood at which the predator feeds on the prey
"""
@inline function predation_loss_preferential(
    P, Z, maximum_grazing_rate, prey_half_saturation, palatability
)
    return maximum_grazing_rate * palatability * holling_type_2(P, prey_half_saturation) * Z
end

"""
    predation_gain_preferential(P, Z, assimilation_efficiency, maximum_grazing_rate, 
        prey_half_saturation, palatability)

Estimates the gain rate of Z (predator) feeding on P (prey).
In this formulation predation rate is multiplied by an assimilation efficiency (β) which
represents 'sloppy feeding'.

!!! formulation

    β * ``g``
    
    where:
    - β = assimilation efficiency of prey to the predator
    - ``g`` = `predation_loss_preferential(P, Z, gₘₐₓ, kₚ, η)`
    - P = phytoplankton concentration
    - Z = zooplankton concentration
    - gₘₐₓ = maximum grazing rate of the predator
    - kₚ = prey half saturation
    - η = palatability

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `assimilation_efficiency`: assimilation efficiency
- `maximum_grazing_rate`: maximum grazing rate
- `prey_half_saturation`: prey density at which predation is half it's maximum rate
- `palatability`: the likelihood at which the predator feeds on the prey
"""
@inline function predation_gain_preferential(
    P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation, palatability
)
    return assimilation_efficiency * predation_loss_preferential(
        P, Z, maximum_grazing_rate, prey_half_saturation, palatability
    )
end

"""
    predation_assimilation_loss_preferential(P, Z, assimilation_efficiency,
        maximum_grazing_rate, prey_half_saturation, palatability)

Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

!!! formulation
    (1-β) * ``g``
    
    where:
    - β = assimilation efficiency of prey to the predator
    - ``g`` = `predation_loss_preferential(P, Z, gₘₐₓ, kₚ, η)`
    - P = phytoplankton concentration
    - Z = zooplankton concentration
    - gₘₐₓ = maximum grazing rate of the predator
    - kₚ = prey half saturation
    - η = palatability

!!! info    
    This function differs from `predation_gain_preferential` as it represents the transfer of 
    biomass from the prey to the environment rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `assimilation_efficiency`: assimilation efficiency
- `maximum_grazing_rate`: maximum grazing rate
- `prey_half_saturation`: prey density at which predation is half it's maximum rate
- `palatability`: the likelihood at which the predator feeds on the prey

"""
@inline function predation_assimilation_loss_preferential(
    P, Z, assimilation_efficiency, maximum_grazing_rate, prey_half_saturation, palatability
)
    return (1 - assimilation_efficiency) * predation_loss_preferential(
        P, Z, maximum_grazing_rate, prey_half_saturation, palatability
    )
end

"""Binary prey data for emergent assimilation efficiency."""
struct AssimilationPreyParameters
    can_be_eaten::Bool
end

"""Predator data for emergent assimilation efficiency."""
struct AssimilationPredatorParameters{FT<:AbstractFloat}
    can_eat::Bool
    assimilation_efficiency::FT
end

"""
    assimilation_efficiency_emergent_binary(prey_data, predator_data)

Determines the assimilation efficiency of a predator consuming prey, based on binary conditions of edibility.
The function evaluates whether the predator can eat the prey and whether the prey can be consumed, and assigns the assimilation efficiency accordingly.

!!! formulation
    β if ``e_{prey}`` = 1 & ``e_{pred}`` = 1

    0 otherwise

    where:
    - β = assimilation efficiency of prey to the predator
    - ``e_{prey}`` = whether prey can be eaten (binary value)
    - ``e_{pred}`` = whether predator can eat (binary value)

# Arguments
- `prey_data`: A dictionary containing prey-specific data:
  - `can_be_eaten`: A binary value (1 or 0) indicating if the prey can be consumed by the predator.
- `predator_data`: A dictionary containing predator-specific data:
  - `can_eat`: A binary value (1 or 0) indicating if the predator can consume prey.
  - `assimilation_efficiency`: The efficiency with which the predator assimilates nutrients from the prey if the conditions are met.

# Returns
- `assimilation_efficiency`:
  - If `can_eat` is 1 and `can_be_eaten` is 1, returns the predator's `assimilation_efficiency`.
  - Otherwise, returns 0.
"""
@inline function assimilation_efficiency_emergent_binary(
    prey::AssimilationPreyParameters, predator::AssimilationPredatorParameters{FT}
) where {FT<:AbstractFloat}
    if predator.can_eat && prey.can_be_eaten
        return predator.assimilation_efficiency
    end
    return zero(FT)
end



"""
    predation_loss_sum(plankton_syms, prey_sym, prey_idx)

Build an allocation-free `Expr` that sums preferential predation loss of `prey_sym` (at index `prey_idx`)
to *all* predator plankton in `plankton_syms`.

The generated expression expects these runtime containers to be in scope:
- `maximum_predation_rate`
- `holling_half_saturation`
- `palatability_matrix`

This is the expression-level counterpart of `predation_loss_preferential`.
"""
function predation_loss_sum(plankton_syms::AbstractVector{Symbol}, prey_sym::Symbol, prey_idx::Int)
    terms = Expr[]
    for (pred_idx, pred_sym) in enumerate(plankton_syms)
        push!(terms, :(predation_loss_preferential(
            $prey_sym,
            $pred_sym,
            maximum_predation_rate[$pred_idx],
            holling_half_saturation[$pred_idx],
            palatability_matrix[$pred_idx, $prey_idx],
        )))
    end
    return sum_expr(terms)
end

"""
    predation_gain_sum(plankton_syms, predator_sym, predator_idx)

Build an allocation-free `Expr` that sums preferential predation gain of `predator_sym` (at index `predator_idx`)
from all prey plankton in `plankton_syms`.

The generated expression expects:
- `assimilation_efficiency_matrix`
- `maximum_predation_rate`
- `holling_half_saturation`
- `palatability_matrix`
"""
function predation_gain_sum(plankton_syms::AbstractVector{Symbol}, predator_sym::Symbol, predator_idx::Int)
    terms = Expr[]
    for (prey_idx, prey_sym) in enumerate(plankton_syms)
        push!(terms, :(predation_gain_preferential(
            $prey_sym,
            $predator_sym,
            assimilation_efficiency_matrix[$predator_idx, $prey_idx],
            maximum_predation_rate[$predator_idx],
            holling_half_saturation[$predator_idx],
            palatability_matrix[$predator_idx, $prey_idx],
        )))
    end
    return sum_expr(terms)
end

"""
    predation_assimilation_loss_sum(plankton_syms)

Build an allocation-free `Expr` that sums assimilation losses from all predator–prey pairs.

The generated expression expects:
- `assimilation_efficiency_matrix`
- `maximum_predation_rate`
- `holling_half_saturation`
- `palatability_matrix`
"""
function predation_assimilation_loss_sum(plankton_syms::AbstractVector{Symbol})
    terms = Expr[]
    for (pred_idx, pred_sym) in enumerate(plankton_syms)
        for (prey_idx, prey_sym) in enumerate(plankton_syms)
            push!(terms, :(predation_assimilation_loss_preferential(
                $prey_sym,
                $pred_sym,
                assimilation_efficiency_matrix[$pred_idx, $prey_idx],
                maximum_predation_rate[$pred_idx],
                holling_half_saturation[$pred_idx],
                palatability_matrix[$pred_idx, $prey_idx],
            )))
        end
    end
    return sum_expr(terms)
end

end # module
