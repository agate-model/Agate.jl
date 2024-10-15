# phytoplankton growth
menden_limitation(R, k) = R / (k + R)
photosynthetic_growth(N, P, μ₀, kₙ) =  μ₀ * menden_limitation(N, kₙ) * P
# zooplankton growth
holling_type_2(R, k) = R / (k + R)
#predation_loss(P, Z, gₘₐₓ, kₚ, palat) = gₘₐₓ * palat * holling_type_2(P, kₚ) * Z
#predation_gain(P, Z, β, gₘₐₓ, kₚ, palat) = predation_rate(P, Z, gₘₐₓ, kₚ, palat) * β
#predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ, palat) = (1 - β) * predation_rate(P, Z, gₘₐₓ, kₚ, palat)
#summed_predation_gain(i, P, β, gₘₐₓ, kₚ, palat) =  sum([predation_gain(P[j], P[i], β[j,i], gₘₐₓ[j], kₚ[j], palat[j,i]) for j in eachindex(P)])
#summed_predation_assimilation_loss(i, P, β, gₘₐₓ, kₚ, palat) =  sum([predation_gain(P[j], P[i], β[j,i], gₘₐₓ[j], kₚ[j], palat[j,i]) for j in eachindex(P)])

function predation_loss(P, Z, maximum_predation_rate, holling_half_saturation, palatability) 
    loss =  maximum_predation_rate * palatability * holling_type_2(P, holling_half_saturation) * Z
    return loss
end

function predation_gain(P, Z, assimilation_efficiency, maximum_predation_rate, holling_half_saturation, palatability) 
    gain = predation_loss(P, Z, maximum_predation_rate, holling_half_saturation, palatability) * assimilation_efficiency
    return gain
end

function predation_assimilation_loss(P, Z, assimilation_efficiency, maximum_predation_rate, holling_half_saturation, palatability) 
    assimilation_loss = predation_loss(P, Z, maximum_predation_rate, holling_half_saturation, palatability) * (1-assimilation_efficiency)
    return assimilation_loss
end

function summed_predation_loss(prey_index, P, maximum_predation_rate, holling_half_saturation, palatability)
    loss = sum(predation_loss(
        P[prey_index], 
        P[predator_index], 
        maximum_predation_rate[predator_index], 
        holling_half_saturation[predator_index], 
        palatability[predator_index, prey_index]
    ) for predator_index in eachindex(P))
    
    return loss
end

function summed_predation_gain(prey_index, P, assimilation_efficiency, maximum_predation_rate, holling_half_saturation, palatability)
    gain = sum(predation_gain(
        P[prey_index], 
        P[predator_index], 
        assimilation_efficiency[predator_index],
        maximum_predation_rate[predator_index], 
        holling_half_saturation[predator_index], 
        palatability[predator_index, prey_index]
    ) for predator_index in eachindex(P))
    
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
net_photosynthetic_growth(N, P, μ₀, kₙ) = sum([photosynthetic_growth(N, P[i], μ₀[i], kₙ[i]) for i in eachindex(P)])
net_predation_assimilation_loss(P, β, gₘₐₓ, kₚ, palat) = sum([summed_predation_assimilation_loss(i, P, β[j,i], gₘₐₓ[j], kₚ[j], palat[j,i]) * β[i] for i in eachindex(P)])
#generic plankton
function plankton_dt(i, N, P, μ₀, kₙ, gₘₐₓ, β, palat)
    growth =
    + photosynthetic_growth(N, P[i], μ₀[i], kₙ[i])
    - plankton_mortality_loss(P[i], l[i])
    - summed_predation_loss(i, P, Z, gₘₐₓ, kₚ, palat)
    + summed_predation_gain(i, P, Z, β, gₘₐₓ, kₚ, palat)
    return growth
end