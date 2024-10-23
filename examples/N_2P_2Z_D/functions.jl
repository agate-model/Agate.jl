# phytoplankton growth
function menden_limitation(N::Float64, nitrogen_half_saturation::Float64)
    return N / (nitrogen_half_saturation + N)
end
function smith_light_limitation(PAR::Float64, alpha::Float64, maximum_growth_rate::Float64)
    return alpha * PAR / sqrt(maximum_growth_rate^2 + alpha^2 * PAR^2)
end
function photosynthetic_growth(
    N::Float64,
    P::Float64,
    PAR::Float64,
    maximum_growth_rate::Float64,
    nitrogen_half_saturation::Float64,
    alpha::Float64,
)
    return maximum_growth_rate *
           menden_limitation(N, nitrogen_half_saturation) *
           smith_light_limitation(PAR, alpha, maximum_growth_rate) *
           P
end
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
    prey_index::Int,
    P::Vector{Float64},
    maximum_predation_rate::Vector{Float64},
    holling_half_saturation::Vector{Float64},
    palatability::Matrix{Float64},
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
    predator_index::Int,
    P::Vector{Float64},
    assimilation_efficiency::Matrix{Float64},
    maximum_predation_rate::Vector{Float64},
    holling_half_saturation::Vector{Float64},
    palatability::Matrix{Float64},
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
estimate the total assimlation loss during predation.

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
- `assimilation_loss`: The summed rate of predation gain for plankton[predator_index]

"""
function summed_predation_assimilation_loss(
    predator_index::Int,
    P::Vector{Float64},
    assimilation_efficiency::Matrix{Float64},
    maximum_predation_rate::Vector{Float64},
    holling_half_saturation::Vector{Float64},
    palatability::Matrix{Float64},
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

#mortality
linear_loss(P::Float64, linear_mortality::Float64) = linear_mortality * P
quadratic_loss(P::Float64, quadratic_mortality::Float64) = quadratic_mortality * P^2
#detritus
function remineralization(D::Float64, detritus_remineralization::Float64)
    return D * detritus_remineralization
end
#sums
"""
Net loss of all plankton due to linear mortality.

# Arguments
- `P::Vector{Symbol}`: Vector which includes all plankton.
- `linear_mortality::Vector{Float}`: Vector of all plankton linear mortality rates.

"""
function net_linear_loss(
    P::Vector{Float64}, linear_mortality::Vector{Float64}, fraction::Float64
)
    return sum([linear_loss(P[i], linear_mortality[i]) for i in eachindex(P)]) * fraction
end
"""
Net loss of all plankton due to quadratic mortality.

# Arguments
- `P::Vector{Symbol}`: Vector which includes all plankton.
- `linear_mortality::Vector{Float}`: Vector of all plankton quadratic mortality rates.

"""
function net_quadratic_loss(
    P::Vector{Float64}, quadratic_mortality::Vector{Float64}, fraction::Float64
)
    return sum(
        [quadratic_loss(P[i], quadratic_mortality[i]) for i in eachindex(P)] * fraction
    )
end
"""
Net photosynthetic growth of all plankton.

# Arguments
- `N::Symbol`: Nitrogen
- `P::Vector{Symbol}`: Vector which includes all plankton.
- `PAR::Symbol`: PAR
- `maximum_growth_rate::Vector{Float}`: Vector of all plankton maximum growth rates.
- `nitrogen_half_saturation::Vector{Float}`: Vector of all plankton nitrogen half saturation constants.

"""
function net_photosynthetic_growth(
    N::Float64,
    P::Vector{Float64},
    PAR::Float64,
    maximum_growth_rate::Vector{Float64},
    nitrogen_half_saturation::Vector{Float64},
    alpha::Vector{Float64},
)
    return sum([
        photosynthetic_growth(
            N, P[i], PAR, maximum_growth_rate[i], nitrogen_half_saturation[i], alpha[i]
        ) for i in eachindex(P)
    ])
end

"""
Net predator assimilation loss of all plankton.

# Arguments
- `P::Vector{Symbol}`: Vector which includes all plankton.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton maximum growth rates.
- `maximum_predation_rate::Vector{Float}`: Vector of all plankton maximum predation rates.
- `palatability::Matrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]). 
    For a non-predator [i,:]=0. 
- `assimilation efficiency::Matrix{Float}`: Matrix of all plankton assimilation efficiencies where:
    each row is a predator and each column is a prey (palat[predator, prey]). 
    For a non-predator [i,:]=0. 
"""
function net_predation_assimilation_loss(
    P::Vector{Float64},
    holling_half_saturation::Vector{Float64},
    maximum_predation_rate::Vector{Float64},
    assimilation_efficiency::Matrix{Float64},
    palatability::Matrix{Float64},
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
- `N::Symbol`: Nitrogen
- `P::Vector{Symbol}`: Vector which includes all plankton.
- `linear_mortality::Vector{Float}`: Vector of all plankton linear mortality rates.
- `maximum_growth_rate::Vector{Float}`: Vector of all plankton maximum growth rates.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton maximum growth rates.
- `maximum_predation_rate::Vector{Float}`: Vector of all plankton maximum predation rates.
- `assimilation efficiency::Matrix{Float}`: Matrix of all plankton assimilation efficiencies where:
    each row is a predator and each column is a prey (palat[predator, prey]). 
    For a non-predator [i,:]=0. 
- `palatability::Matrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]). 
    For a non-predator [i,:]=0. 
"""
function plankton_dt(
    plankton_index::Int,
    N::Float64,
    P::Vector{Float64},
    PAR::Float64,
    linear_mortality::Vector{Float64},
    quadratic_mortality::Vector{Float64},
    maximum_growth_rate::Vector{Float64},
    holling_half_saturation::Vector{Float64},
    nitrogen_half_saturation::Vector{Float64},
    alpha::Vector{Float64},
    maximum_predation_rate::Vector{Float64},
    assimilation_efficiency::Matrix{Float64},
    palatability::Matrix{Float64},
)
    growth =
        photosynthetic_growth(
            N,
            P[plankton_index],
            PAR,
            maximum_growth_rate[plankton_index],
            nitrogen_half_saturation[plankton_index],
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
