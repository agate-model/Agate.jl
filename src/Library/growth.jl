"""
A library of modules to define photosynthetic growth

"""
module Growth

"
    PCⱼ = PCⱼᵐᵃˣ * γⱼⁿᵘᵗ *  γⱼˡⁱᵍʰᵗ * fⱼᵗᵉᵐᵖ *  γⱼᶜᵒ²

Carbon-specific growth rate for plankton j (Default MITgcm-DARWIN formulation).

Where: 
PCᵐᵃˣⱼ = maximum carbon-specific growth rate for plankton j,
γⁿᵘᵗⱼ = nutrient limition,
γˡⁱᵍʰᵗⱼ = light limition,
fᵗᵉᵐᵖⱼ = temperature limitation,
γᶜᵒ²ⱼ = carbon dioxide limitation.

#Is this a useful function?
#could set default to 1 for all limitations?

"
function default_carbon_growth(PCⱼᵐᵃˣ, γⱼⁿᵘᵗ,  γⱼˡⁱᵍʰᵗ, fⱼᵗᵉᵐᵖ,  γⱼᶜᵒ²)
    PCⱼᵐᵃˣ * γⱼⁿᵘᵗ *  γⱼˡⁱᵍʰᵗ * fⱼᵗᵉᵐᵖ *  γⱼᶜᵒ²
end

"
    if Iₜₒₜ > Iₘᵢₙ
        PCᵐⱼ*(1-exp(-(γᶜᶠᵉⱼ*αⱼ*I*Chl¨Cⱼ)/(PCᵐⱼ)))*γⁱⁿʰⱼ
    else
        PCⱼ = 0

Carbon specific growth of plankton j (Geider MITgcm-DARWIN formulation). 

Where: 
Iₜₒₜ = total irradiance accross wavelenghts
Iₘᵢₙ = minimum irradiance for photosynthesis to occur 
γᶜᶠᵉⱼ = effect of iron limitation on photosynthesis of plankton j
Chl¨Cⱼ = chlorophyll to carbon ratio of plankton j
γⁱⁿʰⱼ = light inhibition of plankton j
fᵃᴵⱼ = summed response to irradiance for all wavelenghts of plankton j

"
function geider_carbon_growth(PCᵐⱼ, γⱼⁱⁿʰ, Iₜₒₜ, Iₘᵢₙ, Chl¨Cⱼ, γⱼᶜᶠᵉ, fⱼᵃᴵ)
    if Iₜₒₜ > Iₘᵢₙ
        return PCᵐⱼ*(1-exp(-(γⱼᶜᶠᵉ*fⱼᵃᴵ*Chl¨Cⱼ)/(PCᵐⱼ)))*γⱼⁱⁿʰ
    else
        return 0
    end
end

"
    PCᵐⱼ = PCⱼᵐᵃˣ *  γⱼⁿᵘᵗ * fⱼᵗᵉᵐᵖ *  γⱼᶜᵒ²

Non-light-limited carbon-specific growth of plankton j (Geider MITgcm-DARWIN formulation). 

Where: 
PCᵐᵃˣⱼ = maximum carbon-specific growth rate for plankton j,
γⁿᵘᵗⱼ = nutrient limition,
fᵗᵉᵐᵖⱼ = temperature limitation,
γᶜᵒ²ⱼ = carbon dioxide limitation.

"
function geider_light_saturated_growth(PCⱼᵐᵃˣ, γⱼⁿᵘᵗ, fⱼᵗᵉᵐᵖ,  γⱼᶜᵒ²)
    PCᵐⱼ = PCⱼᵐᵃˣ *  γⱼⁿᵘᵗ * fⱼᵗᵉᵐᵖ *  γⱼᶜᵒ²
end

export
    default_carbon_growth
    geider_carbon_growth
    geider_light_saturated_growth
end #module