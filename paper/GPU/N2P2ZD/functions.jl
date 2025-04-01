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

# Arguments
- `P::Float64`: Prey plankton.
- `Z::Float64`: Predator plankton.
- `maximum_predation_rate::Float`: maximum predation rate of the predator
- `holling_half_saturation::Float`: holling half saturation constant of the predator.
- `palatability::Float`: palatability of prey to the predator.

# Returns
- `loss`: The rate of predation loss of P to Z.
"""
function predation_loss(
    P::Float64,
    Z::Float64,
    maximum_predation_rate::Float64,
    holling_half_saturation::Float64,
    palatability::Float64,
)
    loss =
        maximum_predation_rate *
        palatability *
        holling_type_2(P, holling_half_saturation) *
        Z
    return loss
end

"""
Estimates the gain rate of Z (predator) feeding on P (prey).

# Arguments
- `P::Float64`: Prey plankton.
- `Z::Float64`: Predator plankton.
- `maximum_predation_rate::Float`: maximum predation rate of the predator
- `holling_half_saturation::Float`: holling half saturation constant of the predator.
- `palatability::Float`: palatability of prey to the predator.
- `assimilation_efficiency::Float`: assimilation efficiency of prey to the predator.

# Returns
- `gain`: The rate of predation gain of Z to P.
"""
function predation_gain(
    P::Float64,
    Z::Float64,
    assimilation_efficiency::Float64,
    maximum_predation_rate::Float64,
    holling_half_saturation::Float64,
    palatability::Float64,
)
    gain =
        predation_loss(
            P, Z, maximum_predation_rate, holling_half_saturation, palatability
        ) * assimilation_efficiency
    return gain
end

"""
Estimates the rate at which plankton predation gain is lost due to inefficient assimilation
efficiency (e.g. 'sloppy feeding').

Usually, this is used to estimate fluxes from predation to dissolved and particulate
organic matter (DOM and POM).

# Arguments
- `P::Float64`: Prey plankton.
- `Z::Float64`: Predator plankton.
- `maximum_predation_rate::Float`: maximum predation rate of the predator
- `holling_half_saturation::Float`: holling half saturation constant of the predator.
- `palatability::Float`: palatability of prey to the predator.
- `assimilation_efficiency::Float`: assimilation efficiency of prey to the predator.

# Returns
- `assimilation_loss`: The rate at which predation gain is lost to the environment.
"""
function predation_assimilation_loss(
    P::Float64,
    Z::Float64,
    assimilation_efficiency::Float64,
    maximum_predation_rate::Float64,
    holling_half_saturation::Float64,
    palatability::Float64,
)
    assimilation_loss =
        predation_loss(
            P, Z, maximum_predation_rate, holling_half_saturation, palatability
        ) * (1 - assimilation_efficiency)
    return assimilation_loss
end

"""
Estimates the total loss rate of the prey (P[prey_index]) to predation.

For plankton P[prey_index], the function loops over each predator to estimate the total loss
of plankton i due to predation.

# Arguments
- `prey_index::Int`: Index of the prey, e.g. P[prey_index].
- `P::AbstractVector{Float64}`: AbstractVector which includes all plankton.
- `maximum_predation_rate::AbstractVector{Float}`: AbstractVector of all plankton predation rates.
- `holling_half_saturation::AbstractVector{Float}`: AbstractVector of all plankton predation half saturation constants.
- `palatability::AbstractMatrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.

# Returns
- `loss`: The summed rate of predation loss for plankton[prey_index]
"""
function summed_predation_loss(
    prey_index::Int,
    P::AbstractVector{Float64},
    maximum_predation_rate::AbstractVector{Float64},
    holling_half_saturation::AbstractVector{Float64},
    palatability::AbstractMatrix{Float64},
)
    loss = sum(
        predation_loss(
            P[prey_index],
            P[predator_index],
            maximum_predation_rate[predator_index],
            holling_half_saturation[predator_index],
            palatability[predator_index, prey_index],
        ) for predator_index in eachindex(P)
    )

    return loss
end

"""
Estimates the total predation gain of the predator (P[predator_index]) feeding on all plankton.

For plankton P[predator_index], the function loops over each prey (P[prey_index]) to
estimate the total gain due to predation.

# Arguments
- `predator_index::Int`: Index of the predator, e.g. P[predator_index].
- `P::AbstractVector{Float64}`: AbstractVector which includes all plankton.
- `maximum_predation_rate::AbstractVector{Float}`: AbstractVector of all plankton predation rates.
- `holling_half_saturation::AbstractVector{Float}`: AbstractVector of all plankton predation half saturation constants.
- `palatability::AbstractMatrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.
- `assimilation efficiency::AbstractMatrix{Float}`: Matrix of all plankton assimilation efficiencies where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.

# Returns
- `gain`: The summed rate of predation gain for plankton[predator_index]
"""
function summed_predation_gain(
    predator_index::Int,
    P::AbstractVector{Float64},
    assimilation_efficiency::AbstractMatrix{Float64},
    maximum_predation_rate::AbstractVector{Float64},
    holling_half_saturation::AbstractVector{Float64},
    palatability::AbstractMatrix{Float64},
)
    gain = sum(
        predation_gain(
            P[prey_index],
            P[predator_index],
            assimilation_efficiency[predator_index, prey_index],
            maximum_predation_rate[predator_index],
            holling_half_saturation[predator_index],
            palatability[predator_index, prey_index],
        ) for prey_index in eachindex(P)
    )

    return gain
end

"""
Estimates the total assimilation loss of the predator (P[predator_index]) feeding on all plankton.

For plankton P[predator_index], the function loops over each prey (P[prey_index]) to
estimate the total assimilation loss during predation.

# Arguments
- `predator_index::Int`: Index of the predator, e.g. P[predator_index].
- `P::AbstractVector{Float64}`: AbstractVector which includes all plankton.
- `maximum_predation_rate::AbstractVector{Float}`: AbstractVector of all plankton predation rates.
- `holling_half_saturation::AbstractVector{Float}`: AbstractVector of all plankton predation half saturation constants.
- `palatability::AbstractMatrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.
- `assimilation efficiency::AbstractMatrix{Float}`: Matrix of all plankton assimilation efficiencies where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.

# Returns
- `assimilation_loss`: The summed rate of predation gain for plankton[predator_index]
"""
function summed_predation_assimilation_loss(
    predator_index::Int,
    P::AbstractVector{Float64},
    assimilation_efficiency::AbstractMatrix{Float64},
    maximum_predation_rate::AbstractVector{Float64},
    holling_half_saturation::AbstractVector{Float64},
    palatability::AbstractMatrix{Float64},
)
    assimilation_loss = sum(
        predation_assimilation_loss(
            P[prey_index],
            P[predator_index],
            assimilation_efficiency[predator_index, prey_index],
            maximum_predation_rate[predator_index],
            holling_half_saturation[predator_index],
            palatability[predator_index, prey_index],
        ) for prey_index in eachindex(P)
    )

    return assimilation_loss
end

#sums
"""
 Net loss of all plankton due to linear mortality.

 # Arguments
 - `P::AbstractVector{Float64}`: AbstractVector which includes all plankton.
 - `linear_mortality::AbstractVector{Float}`: AbstractVector of all plankton linear mortality rates.
 """
function custom_net_linear_loss(
    P::AbstractVector{Float64}, linear_mortality::AbstractVector{Float64}, fraction::Float64
)
    return sum([linear_loss(P[i], linear_mortality[i]) for i in eachindex(P)]) * fraction
end

"""
Net loss of all plankton due to quadratic mortality.

# Arguments
- `P::Vector{Float64}`: Vector which includes all plankton.
- `linear_mortality::Vector{Float}`: Vector of all plankton quadratic mortality rates.
"""
function custom_net_quadratic_loss(
    P::AbstractVector{Float64}, quadratic_mortality::AbstractVector{Float64}, fraction::Float64
)
    return sum(
        [quadratic_loss(P[i], quadratic_mortality[i]) for i in eachindex(P)] * fraction
    )
end

"""
Net photosynthetic growth of all plankton.

# Arguments
- `N::Float64`: nutrient
- `P::Vector{Float64}`: Vector which includes all plankton.
- `PAR::Float64`: PAR
- `maximum_growth_rate::Vector{Float}`: Vector of all plankton maximum growth rates.
- `nutrient_half_saturation::Vector{Float}`: Vector of all plankton nutrient half saturation constants.
"""
function net_photosynthetic_growth(
    N::Float64,
    P::AbstractVector{Float64},
    PAR::Float64,
    maximum_growth_rate::AbstractVector{Float64},
    nutrient_half_saturation::AbstractVector{Float64},
    alpha::AbstractVector{Float64},
)
    return sum([
        photosynthetic_growth(
            N, P[i], PAR, maximum_growth_rate[i], nutrient_half_saturation[i], alpha[i]
        ) for i in eachindex(P)
    ])
end

"""
Net predator assimilation loss of all plankton.

# Arguments
- `P::Vector{Float64}`: Vector which includes all plankton.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton predation half saturation constants.
- `maximum_predation_rate::Vector{Float}`: Vector of all plankton maximum predation rates.
- `palatability::AbstractMatrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.
- `assimilation efficiency::AbstractMatrix{Float}`: Matrix of all plankton assimilation efficiencies where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.
"""
function net_predation_assimilation_loss(
    P::AbstractVector{Float64},
    holling_half_saturation::AbstractVector{Float64},
    maximum_predation_rate::AbstractVector{Float64},
    assimilation_efficiency::AbstractMatrix{Float64},
    palatability::AbstractMatrix{Float64},
)
    return sum([
        summed_predation_assimilation_loss(
            predator_index,
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
        ) for predator_index in eachindex(P)
    ])
end
#generic plankton

"""
Wrapper function to estimate the rate at which plankton biomass changes over time.

# Arguments
- `plankton_index:Int`: The index of the plankton for which the rate of change is estimated
- `N::Float64`: nutrient
- `P::AbstractVector{Float64}`: AbstractVector which includes all plankton.
- `linear_mortality::AbstractVector{Float}`: AbstractVector of all plankton linear mortality rates.
- `maximum_growth_rate::AbstractVector{Float}`: AbstractVector of all plankton maximum growth rates.
- `holling_half_saturation::AbstractVector{Float}`: AbstractVector of all plankton predation half saturation constants.
- `maximum_predation_rate::AbstractVector{Float}`: AbstractVector of all plankton maximum predation rates.
- `assimilation efficiency::AbstractMatrix{Float}`: Matrix of all plankton assimilation efficiencies where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.
- `palatability::AbstractMatrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.
"""
function plankton_dt(
    plankton_index::Int,
    N::Float64,
    P::AbstractVector{Float64},
    PAR::Float64,
    linear_mortality::AbstractVector{Float64},
    quadratic_mortality::AbstractVector{Float64},
    maximum_growth_rate::AbstractVector{Float64},
    holling_half_saturation::AbstractVector{Float64},
    nutrient_half_saturation::AbstractVector{Float64},
    alpha::AbstractVector{Float64},
    maximum_predation_rate::AbstractVector{Float64},
    assimilation_efficiency::AbstractMatrix{Float64},
    palatability::AbstractMatrix{Float64},
)
    growth =
        photosynthetic_growth(
            N,
            P[plankton_index],
            PAR,
            maximum_growth_rate[plankton_index],
            nutrient_half_saturation[plankton_index],
            alpha[plankton_index],
        ) - linear_loss(P[plankton_index], linear_mortality[plankton_index]) -
        quadratic_loss(P[plankton_index], quadratic_mortality[plankton_index]) -
        summed_predation_loss(
            plankton_index, P, maximum_predation_rate, holling_half_saturation, palatability
        ) + summed_predation_gain(
            plankton_index,
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
        )
    return growth
end
