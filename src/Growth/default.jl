"
    PCⱼ = PCᵐᵃˣⱼ * γⁿᵘᵗⱼ *  γˡⁱᵍʰᵗⱼ * fᵗᵉᵐᵖⱼ *  γᶜᵒ²ⱼ

Carbon-specific growth rate for plankton j (Default MITgcm-DARWIN formulation).

Where: 
PCᵐᵃˣⱼ = maximum carbon-specific growth rate for plankton j
γⁿᵘᵗⱼ = nutrient limition
γˡⁱᵍʰᵗⱼ = light limition
fᵗᵉᵐᵖⱼ = temperature limitation
γᶜᵒ²ⱼ = carbon dioxide limitation

"
function default_PCⱼ(PCᵐᵃˣⱼ, γⁿᵘᵗⱼ,  γˡⁱᵍʰᵗⱼ, fᵗᵉᵐᵖⱼ,  γᶜᵒ²ⱼ)
    PCᵐᵃˣⱼ * γⁿᵘᵗⱼ *  γˡⁱᵍʰᵗⱼ * fᵗᵉᵐᵖⱼ *  γᶜᵒ²ⱼ
end


"
    γˡⁱᵍʰᵗⱼ = (1 - ℯ^(kˢᵃᵗⱼ*I)) * ℯ^kⁱⁿʰⱼ * nˡⁱᵍʰᵗⱼ

Light limitation for plankton j (Default MITgcm-DARWIN formulation). 

Where: 
kˢᵃᵗⱼ = half saturation constant of light saturation of plankton j.
I = irradiance
kⁱⁿʰⱼ = half saturation constant of light inhibition of plankton j.
nˡⁱᵍʰᵗⱼ = light penalty term of plankton j

"
function γˡⁱᵍʰᵗⱼ(I, kˢᵃᵗⱼ, kⁱⁿʰⱼ, nˡⁱᵍʰᵗⱼ)
    (1 - ℯ^(kˢᵃᵗⱼ*I)) * ℯ^kⁱⁿʰⱼ * nˡⁱᵍʰᵗⱼ
end