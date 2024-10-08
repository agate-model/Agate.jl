"""
Modules related to phytoplankton light uptake

"""
module Photosynthesis

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
    return (1 - ℯ^(kⱼˢᵃᵗ * I)) * ℯ^kⱼⁱⁿʰ * nⱼˡⁱᵍʰᵗ
end

"
    α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)

Smith 1936 formulation of light limitation (also see Evans and Parslow, 1985).

Where:
α = Initial photosynthetic slope
PAR = Photosynthetic Active Radiation
μ₀ = Maximum growth rate at T = 0 °C (this seems weird?, from Kuhn 2015)

"
function smith_light_limitation(PAR, α, μ₀)
    return α * PAR / sqrt(μ₀^2 + α^2 * PAR^2)
end

export γⱼˡⁱᵍʰᵗ
smith_light_limitation
end # module
