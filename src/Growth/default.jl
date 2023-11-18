"
    PCⱼ = PCⱼᵐᵃˣ * γⱼⁿᵘᵗ *  γⱼˡⁱᵍʰᵗ * fⱼᵗᵉᵐᵖ *  γⱼᶜᵒ²

Carbon-specific growth rate for plankton j.

Where: 
PCᵐᵃˣⱼ = maximum carbon-specific growth rate for plankton j
γⁿᵘᵗⱼ = nutrient limition
γˡⁱᵍʰᵗⱼ = light limition
fᵗᵉᵐᵖⱼ = temperature limitation
γᶜᵒ²ⱼ = carbon dioxide limitation

"
function PCⱼ(PCⱼᵐᵃˣ, γⱼⁿᵘᵗ,  γⱼˡⁱᵍʰᵗ, fⱼᵗᵉᵐᵖ,  γⱼᶜᵒ²)
    PCⱼᵐᵃˣ * γⱼⁿᵘᵗ *  γⱼˡⁱᵍʰᵗ * fⱼᵗᵉᵐᵖ *  γⱼᶜᵒ²
end


"
    γⱼˡⁱᵍʰᵗ = (1 - ℯ^(kⱼˢᵃᵗ*I)) * ℯ^kⱼⁱⁿʰ * nⱼˡⁱᵍʰᵗ

Light limitation for plankton j (Default MITgcm-DARWIN formulation). 

Where: 
kˢᵃᵗⱼ = half saturation constant of light saturation of plankton j.
I = irradiance
kⁱⁿʰⱼ = half saturation constant of light inhibition of plankton j.
nˡⁱᵍʰᵗⱼ = light penalty term of plankton j

"
function default_γⱼˡⁱᵍʰᵗ(I, kⱼˢᵃᵗ, kⱼⁱⁿʰ, nⱼˡⁱᵍʰᵗ)
    (1 - ℯ^(kⱼˢᵃᵗ*I)) * ℯ^kⱼⁱⁿʰ * nⱼˡⁱᵍʰᵗ
end