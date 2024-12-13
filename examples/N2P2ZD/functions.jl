function remineralization(D::Real, detritus_remineralization::Real)
    return D * detritus_remineralization
end

# phytoplankton growth
function photosynthetic_growth(
    N::Real,
    P::Real,
    PAR::Real,
    maximum_growth_rate::Real,
    nitrogen_half_saturation::Real,
    alpha::Real,
)
    return maximum_growth_rate *
           monod_limitation(N, nitrogen_half_saturation) *
           smith_light_limitation(PAR, alpha, maximum_growth_rate) *
           P
end
# zooplankton growth
holling_type_2(R::Real, k::Real) = R / (k + R)

"""
Estimates the loss rate of P (prey), to Z (predator).

# Arguments
- `P::Real`: Prey plankton.
- `Z::Real`: Predator plankton.
- `maximum_predation_rate::Float`: maximum predation rate of the predator
- `holling_half_saturation::Float`: holling half saturation constant of the predator.
- `palatability::Float`: palatability of prey to the predator.

# Returns
- `loss`: The rate of predation loss of P to Z.
"""
function predation_loss(
    P::Real,
    Z::Real,
    maximum_predation_rate::Real,
    holling_half_saturation::Real,
    palatability::Real,
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
- `P::Real`: Prey plankton.
- `Z::Real`: Predator plankton.
- `maximum_predation_rate::Float`: maximum predation rate of the predator
- `holling_half_saturation::Float`: holling half saturation constant of the predator.
- `palatability::Float`: palatability of prey to the predator.
- `assimilation_efficiency::Float`: assimilation efficiency of prey to the predator.

# Returns
- `gain`: The rate of predation gain of Z to P.
"""
function predation_gain(
    P::Real,
    Z::Real,
    assimilation_efficiency::Real,
    maximum_predation_rate::Real,
    holling_half_saturation::Real,
    palatability::Real,
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
- `P::Real`: Prey plankton.
- `Z::Real`: Predator plankton.
- `maximum_predation_rate::Float`: maximum predation rate of the predator
- `holling_half_saturation::Float`: holling half saturation constant of the predator.
- `palatability::Float`: palatability of prey to the predator.
- `assimilation_efficiency::Float`: assimilation efficiency of prey to the predator.

# Returns
- `assimilation_loss`: The rate at which predation gain is lost to the environment.
"""
function predation_assimilation_loss(
    P::Real,
    Z::Real,
    assimilation_efficiency::Real,
    maximum_predation_rate::Real,
    holling_half_saturation::Real,
    palatability::Real,
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
- `P::Vector{<:Real}`: Vector which includes all plankton.
- `maximum_predation_rate::Vector{Float}`: Vector of all plankton predation rates.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton predation half saturation constants.
- `palatability::Matrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.

# Returns
- `loss`: The summed rate of predation loss for plankton[prey_index]
"""
function summed_predation_loss(
    prey_index::Int,
    P::Vector{<:Real},
    maximum_predation_rate::Vector{<:Real},
    holling_half_saturation::Vector{<:Real},
    palatability::Matrix{<:Real},
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
- `P::Vector{<:Real}`: Vector which includes all plankton.
- `maximum_predation_rate::Vector{Float}`: Vector of all plankton predation rates.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton predation half saturation constants.
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
    P::Vector{<:Real},
    assimilation_efficiency::Matrix{<:Real},
    maximum_predation_rate::Vector{<:Real},
    holling_half_saturation::Vector{<:Real},
    palatability::Matrix{<:Real},
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
- `P::Vector{<:Real}`: Vector which includes all plankton.
- `maximum_predation_rate::Vector{Float}`: Vector of all plankton predation rates.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton predation half saturation constants.
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
    P::Vector{<:Real},
    assimilation_efficiency::Matrix{<:Real},
    maximum_predation_rate::Vector{<:Real},
    holling_half_saturation::Vector{<:Real},
    palatability::Matrix{<:Real},
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
 - `P::Vector{<:Real}`: Vector which includes all plankton.
 - `linear_mortality::Vector{Float}`: Vector of all plankton linear mortality rates.

 """
function custom_net_linear_loss(
    P::Vector{<:Real}, linear_mortality::Vector{<:Real}, fraction::Real
)
    return sum([linear_loss(P[i], linear_mortality[i]) for i in eachindex(P)]) * fraction
end
"""
Net loss of all plankton due to quadratic mortality.
# Arguments
- `P::Vector{<:Real}`: Vector which includes all plankton.
- `linear_mortality::Vector{Float}`: Vector of all plankton quadratic mortality rates.
"""
function custom_net_quadratic_loss(
    P::Vector{<:Real}, quadratic_mortality::Vector{<:Real}, fraction::Real
)
    return sum(
        [quadratic_loss(P[i], quadratic_mortality[i]) for i in eachindex(P)] * fraction
    )
end

"""
Net photosynthetic growth of all plankton.

# Arguments
- `N::Real`: Nitrogen
- `P::Vector{<:Real}`: Vector which includes all plankton.
- `PAR::Real`: PAR
- `maximum_growth_rate::Vector{Float}`: Vector of all plankton maximum growth rates.
- `nitrogen_half_saturation::Vector{Float}`: Vector of all plankton nitrogen half saturation constants.
"""
function net_photosynthetic_growth(
    N::Real,
    P::Vector{<:Real},
    PAR::Real,
    maximum_growth_rate::Vector{<:Real},
    nitrogen_half_saturation::Vector{<:Real},
    alpha::Vector{<:Real},
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
- `P::Vector{<:Real}`: Vector which includes all plankton.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton predation half saturation constants.
- `maximum_predation_rate::Vector{Float}`: Vector of all plankton maximum predation rates.
- `palatability::Matrix{Float}`: Matrix of all plankton palatabilities where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.
- `assimilation efficiency::Matrix{Float}`: Matrix of all plankton assimilation efficiencies where:
    each row is a predator and each column is a prey (palat[predator, prey]).
    For a non-predator [i,:]=0.
"""
function net_predation_assimilation_loss(
    P::Vector{<:Real},
    holling_half_saturation::Vector{<:Real},
    maximum_predation_rate::Vector{<:Real},
    assimilation_efficiency::Matrix{<:Real},
    palatability::Matrix{<:Real},
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
- `N::Real`: Nitrogen
- `P::Vector{<:Real}`: Vector which includes all plankton.
- `linear_mortality::Vector{Float}`: Vector of all plankton linear mortality rates.
- `maximum_growth_rate::Vector{Float}`: Vector of all plankton maximum growth rates.
- `holling_half_saturation::Vector{Float}`: Vector of all plankton predation half saturation constants.
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
    N::Real,
    P::Vector{<:Real},
    PAR::Real,
    linear_mortality::Vector{<:Real},
    quadratic_mortality::Vector{<:Real},
    maximum_growth_rate::Vector{<:Real},
    holling_half_saturation::Vector{<:Real},
    nitrogen_half_saturation::Vector{<:Real},
    alpha::Vector{<:Real},
    maximum_predation_rate::Vector{<:Real},
    assimilation_efficiency::Matrix{<:Real},
    palatability::Matrix{<:Real},
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
