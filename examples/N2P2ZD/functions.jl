# phytoplankton growth
menden_limitation(R, k) = R / (k + R)
photosynthetic_growth(N, P, μ₀, kₙ) = μ₀ * menden_limitation(N, kₙ) * P
# zooplankton growth
holling_type_2(R, k) = R / (k + R)

"""
Estimates the loss rate of P (prey), to Z (predator).

# Arguments
- `P::Symbol`: Prey plankton.
- `Z::Symbol`: Predator plankton.
- `maximum_predation_rate::Float`: maximum predation rate of the predator
- `holling_half_saturation::Float`: holling half saturation constant of the predator.
- `palatability::Float`: palatability of prey to the predator.

# Returns
- `loss`: The rate of predation loss of P to Z.

"""
function predation_loss(P, Z, maximum_predation_rate, holling_half_saturation, palatability)
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
- `P::Symbol`: Prey plankton.
- `Z::Symbol`: Predator plankton.
- `maximum_predation_rate::Float`: maximum predation rate of the predator
- `holling_half_saturation::Float`: holling half saturation constant of the predator.
- `palatability::Float`: palatability of prey to the predator.
- `assimilation_efficiency::Float`: assimilation efficiency of prey to the predator.

# Returns
- `gain`: The rate of predation gain of Z to P.

"""
function predation_gain(
    P,
    Z,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    gain =
        predation_loss(
            P, Z, maximum_predation_rate, holling_half_saturation, palatability
        ) * assimilation_efficiency
    return gain
end

"""
Estimates the rate at which plankton predation gain is lost due to inefficient assimilation efficiency 
(e.g. 'sloppy feeding').

Usually, this is used to estimate fluxes from predation to dissolved and particulate 
organic matter (DOM and POM).

# Arguments
- `P::Symbol`: Prey plankton.
- `Z::Symbol`: Predator plankton.
- `maximum_predation_rate::Float`: maximum predation rate of the predator
- `holling_half_saturation::Float`: holling half saturation constant of the predator.
- `palatability::Float`: palatability of prey to the predator.
- `assimilation_efficiency::Float`: assimilation efficiency of prey to the predator.

# Returns
- `assimilation_loss`: The rate at which predation gain is lost to the environment.

"""
function predation_assimilation_loss(
    P,
    Z,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    assimilation_loss =
        predation_loss(
            P, Z, maximum_predation_rate, holling_half_saturation, palatability
        ) * (1 - assimilation_efficiency)
    return assimilation_loss
end

"""
Estimates the total loss rate of the prey (P[prey_index]) to predation.

For plankton P[prey_index], the function loops over each predator to 
estimate the total loss of plankton i due to predation.

# Arguments
- `prey_index::Int`: Index of the prey, e.g. P[prey_index].
- `P::Vector{Symbol}`: Vector which includes all plankton.
- `maximum_predation_rate::Vector{Float}`: Vector of all plankton predation rates.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton predation rates.
- `palatability::Matrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]). 
    For a non-predator [i,:]=0. 

# Returns
- `loss`: The summed rate of predation loss for plankton[prey_index]

"""
function summed_predation_loss(
    prey_index, P, maximum_predation_rate, holling_half_saturation, palatability
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
- `P::Vector{Symbol}`: Vector which includes all plankton.
- `maximum_predation_rate::Vector{Float}`: Vector of all plankton predation rates.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton predation rates.
- `palatability::Matrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]). 
    For a non-predator [i,:]=0. 
- `assimilation efficiency::Matrix{Float}`: Matrix of all plankton assimilation efficiencies where:
    each row is a predator and each column is a prey (palat[predator, prey]). 
    For a non-predator [i,:]=0. 

# Returns
- `gain`: The summed rate of predation gain for plankton[predator_index]

"""
function summed_predation_gain(
    predator_index,
    P,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    gain = sum(
        predation_gain(
            P[prey_index],
            P[predator_index],
            assimilation_efficiency[predator_index],
            maximum_predation_rate[predator_index],
            holling_half_saturation[predator_index],
            palatability[predator_index, prey_index],
        ) for prey_index in eachindex(P)
    )

    return gain
end

#mortality
linear_loss(P, l) = l * P
quadratic_loss(P, l) = l * P^2
#detritus
remineralization(D, detritus_remineralization) = D * detritus_remineralization
#sums
net_linear_loss(P, l) = sum([linear_loss(P[i], l[i]) for i in eachindex(P)])
net_quadratic_loss(P, l) = sum([quadratic_loss(P[i], l[i]) for i in eachindex(P)])
function net_photosynthetic_growth(N, P, μ₀, kₙ)
    return sum([photosynthetic_growth(N, P[i], μ₀[i], kₙ[i]) for i in eachindex(P)])
end
function net_predation_assimilation_loss(P, β, gₘₐₓ, kₚ, palat)
    return sum([
        summed_predation_assimilation_loss(i, P, β[j, i], gₘₐₓ[j], kₚ[j], palat[j, i]) *
        β[i] for i in eachindex(P)
    ])
end
#generic plankton

function plankton_dt(
    plankton_index,
    N,
    P,
    linear_mortality,
    maximum_growth_rate,
    holling_half_saturation,
    maximum_predation_rate,
    assimilation_efficiency,
    palatability,
)
    growth =
        +photosynthetic_growth(
            N,
            P[plankton_index],
            maximum_growth_rate[plankton_index],
            holling_half_saturation[plankton_index],
        )
    -plankton_mortality_loss(P[plankton_index], linear_mortality[plankton_index])
    -summed_predation_loss(
        prey_index, P, maximum_predation_rate, holling_half_saturation, palatability
    )
    +summed_predation_gain(
        predator_index,
        P,
        assimilation_efficiency,
        maximum_predation_rate,
        holling_half_saturation,
        palatability,
    )
    return growth
end