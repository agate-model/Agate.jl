using Agate.Library.Mortality
using Agate.Library.Nutrients
using Agate.Library.Photosynthesis

function remineralization(D::Float64, detritus_remineralization::Float64)
    return D * detritus_remineralization
end

# phytoplankton growth
function photosynthetic_growth(
    N::Float64,
    P::Float64,
    PAR::Float64,
    maximum_growth_rate::Float64,
    nutrient_half_saturation::Float64,
    alpha::Float64,
)
    return maximum_growth_rate *
           monod_limitation(N, nutrient_half_saturation) *
           light_limitation_smith(PAR, alpha, maximum_growth_rate) *
           P
end

# zooplankton growth
holling_type_2(R::Float64, k::Float64) = R / (k + R)

"""
Estimates the loss rate of P (prey), to Z (predator).
"""
function predation_loss(
    P::Float64,
    Z::Float64,
    maximum_predation_rate::Float64,
    holling_half_saturation::Float64,
    palatability::Float64,
)
    return maximum_predation_rate * palatability * holling_type_2(P, holling_half_saturation) * Z
end

"""
Estimates the gain rate of Z (predator) feeding on P (prey).
"""
function predation_gain(
    P::Float64,
    Z::Float64,
    assimilation_efficiency::Float64,
    maximum_predation_rate::Float64,
    holling_half_saturation::Float64,
    palatability::Float64,
)
    return predation_loss(P, Z, maximum_predation_rate, holling_half_saturation, palatability) * assimilation_efficiency
end

"""
Estimates the rate at which plankton predation gain is lost due to inefficient assimilation.
"""
function predation_assimilation_loss(
    P::Float64,
    Z::Float64,
    assimilation_efficiency::Float64,
    maximum_predation_rate::Float64,
    holling_half_saturation::Float64,
    palatability::Float64,
)
    return predation_loss(P, Z, maximum_predation_rate, holling_half_saturation, palatability) * (1 - assimilation_efficiency)
end

"""
Estimates the total loss rate of the prey (P[prey_index]) to predation.
"""
function summed_predation_loss(
    prey_index::Int,
    P::AbstractVector{Float64},
    maximum_predation_rate::AbstractVector{Float64},
    holling_half_saturation::AbstractVector{Float64},
    palatability::AbstractMatrix{Float64},
)
    P_prey = P[prey_index]
    P_pred = P
    rates = maximum_predation_rate
    half_sats = holling_half_saturation
    palats = palatability[:, prey_index]
    
    losses = rates .* palats .* holling_type_2.(P_prey, half_sats) .* P_pred
    return sum(losses)
end

"""
Estimates the total predation gain of the predator (P[predator_index]) feeding on all plankton.
"""
function summed_predation_gain(
    predator_index::Int,
    P::AbstractVector{Float64},
    assimilation_efficiency::AbstractMatrix{Float64},
    maximum_predation_rate::AbstractVector{Float64},
    holling_half_saturation::AbstractVector{Float64},
    palatability::AbstractMatrix{Float64},
)
    P_pred = P[predator_index]
    P_prey = P
    rates = maximum_predation_rate[predator_index]
    half_sats = holling_half_saturation[predator_index]
    palats = palatability[predator_index, :]
    effs = assimilation_efficiency[predator_index, :]
    
    gains = rates .* palats .* effs .* holling_type_2.(P_prey, half_sats) .* P_pred
    return sum(gains)
end

"""
Estimates the total assimilation loss of the predator (P[predator_index]) feeding on all plankton.
"""
function summed_predation_assimilation_loss(
    predator_index::Int,
    P::AbstractVector{Float64},
    assimilation_efficiency::AbstractMatrix{Float64},
    maximum_predation_rate::AbstractVector{Float64},
    holling_half_saturation::AbstractVector{Float64},
    palatability::AbstractMatrix{Float64},
)
    P_pred = P[predator_index]
    P_prey = P
    rates = maximum_predation_rate[predator_index]
    half_sats = holling_half_saturation[predator_index]
    palats = palatability[predator_index, :]
    effs = assimilation_efficiency[predator_index, :]
    
    losses = rates .* palats .* (1 .- effs) .* holling_type_2.(P_prey, half_sats) .* P_pred
    return sum(losses)
end

"""
Net loss of all plankton due to linear mortality.
"""
function custom_net_linear_loss(
    P::AbstractVector{Float64},
    linear_mortality::AbstractVector{Float64},
    fraction::Float64
)
    return sum(linear_loss.(P, linear_mortality)) * fraction
end

"""
Net loss of all plankton due to quadratic mortality.
"""
function custom_net_quadratic_loss(
    P::AbstractVector{Float64},
    quadratic_mortality::AbstractVector{Float64},
    fraction::Float64
)
    return sum(quadratic_loss.(P, quadratic_mortality)) * fraction
end

"""
Net photosynthetic growth of all plankton.
"""
function net_photosynthetic_growth(
    N::Float64,
    P::AbstractVector{Float64},
    PAR::Float64,
    maximum_growth_rate::AbstractVector{Float64},
    nutrient_half_saturation::AbstractVector{Float64},
    alpha::AbstractVector,
)
    return sum(photosynthetic_growth.(
        N, P, PAR, maximum_growth_rate, nutrient_half_saturation, alpha
    ))
end

"""
Net predator assimilation loss of all plankton.
"""
function net_predation_assimilation_loss(
    P::AbstractVector{Float64},
    holling_half_saturation::AbstractVector{Float64},
    maximum_predation_rate::AbstractVector{Float64},
    assimilation_efficiency::AbstractMatrix{Float64},
    palatability::AbstractMatrix{Float64},
)
    # Calculate the denominator for Holling type II functional response
    denominator = holling_half_saturation .+ sum(palatability .* P', dims=2)
    
    # Calculate the predation rates for all predator-prey pairs
    predation_rates = maximum_predation_rate .* P .* (palatability .* P') ./ denominator
    
    # Calculate assimilation losses
    assimilation_losses = assimilation_efficiency .* predation_rates
    
    # Sum all the losses
    return sum(assimilation_losses)
end