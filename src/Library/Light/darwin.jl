"
    γⱼˡⁱᵍʰᵗ = (1 - ℯ^(kⱼˢᵃᵗ*I)) * ℯ^kⱼⁱⁿʰ * nⱼˡⁱᵍʰᵗ

Light limitation for plankton j (Default MITgcm-DARWIN formulation). 

Where: 
kⱼˢᵃᵗ = half saturation constant of light saturation of plankton j,
I = irradiance,
kⱼⁱⁿʰ = half saturation constant of light inhibition of plankton j,
nⱼˡⁱᵍʰᵗ = light penalty term of plankton j

"
function γⱼˡⁱᵍʰᵗ(I, kⱼˢᵃᵗ, kⱼⁱⁿʰ, nⱼˡⁱᵍʰᵗ)
    (1 - ℯ^(kⱼˢᵃᵗ*I)) * ℯ^kⱼⁱⁿʰ * nⱼˡⁱᵍʰᵗ
end


