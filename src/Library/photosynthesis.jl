"""
Modules related to phytoplankton photosynthetic growth

"""
module Photosynthesis

"
    α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)

Smith 1936 formulation of light limitation (also see Evans and Parslow, 1985).

Where: 
α = Initial photosynthetic slope
PAR = Photosynthetic Active Radiation
μ₀ = Maximum growth rate at T = 0 °C (this seems weird?, from Kuhn 2015)

"
function smith_light_limitation(PAR, α, μ₀)
    α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)
end


"
    γⱼˡⁱᵍʰᵗ = (1 - ℯ^(kⱼˢᵃᵗ*I)) * ℯ^kⱼⁱⁿʰ * nⱼˡⁱᵍʰᵗ

Light limitation for plankton j (Default MITgcm-DARWIN formulation). 
    
Where: 
kⱼˢᵃᵗ = half saturation constant of light saturation of plankton j,
I = irradiance,
kⱼⁱⁿʰ = half saturation constant of light inhibition of plankton j,
nⱼˡⁱᵍʰᵗ = light penalty term of plankton j

"
function darwin_default_light_limitation(I, kⱼˢᵃᵗ, kⱼⁱⁿʰ, nⱼˡⁱᵍʰᵗ)
    (1 - ℯ^(kⱼˢᵃᵗ*I)) * ℯ^kⱼⁱⁿʰ * nⱼˡⁱᵍʰᵗ
end

"
    Chl¨Cᵃᶜˡⱼ = (Chl¨Cᵐᵃˣⱼ/((1 + Chl¨Cᵐᵃˣⱼ*aⱼ*I)/(2*PCᵐⱼ)))

Acclimated chlorophyll quota for plankton j 
(computed indepently of chlorophyll and carbon tracers).

Where: 
Chl¨Cᵐⁱⁿⱼ = min Chl to Carbon ratio
Chl¨Cᵐᵃˣⱼ = max Chl to Carbon ratio
I = irradiance
aⱼ = slope of irradiance curve
PCᵐⱼ = realised maximum growth rate (ignoring light)

"
function acclimated_chl_carbon_ratio(Chl¨Cⱼᵐⁱⁿ, Chl¨Cⱼᵐᵃˣ, aⱼ, I, PCⱼᵐ)
    Chl¨Cⱼᵃᶜˡ = (Chl¨Cⱼᵐᵃˣ/((1 + Chl¨Cⱼᵐᵃˣ*aⱼ*I)/(2*PCⱼᵐ)))
    if PCⱼᵐ == 0
        Chl¨Cⱼᵐⁱⁿ = 0
    end
    if Chl¨Cⱼᵃᶜˡ < Chl¨Cⱼᵐⁱⁿ
        return Chl¨Cⱼᵐⁱⁿ
    elseif Chl¨Cⱼᵃᶜˡ > Chl¨Cⱼᵐᵃˣ
        return Chl¨Cⱼᵐᵃˣ
    else
        return Chl¨Cⱼᵃᶜˡ
    end
end

"
    a̅ⱼ = sum(Δλαⱼₗᶜʰˡ)/sum(Δλ)

Summed slope of irradiance curve (across wavelenghts) for plankton j.

Where: 
Δ = 
λ = 
αᶜʰˡⱼₗ = chlorophyll specific slope of irradiance curve

"
function summed_irradiance_curve(Δλαⱼₗᶜʰˡ, Δλ)
    sum(Δλαⱼₗᶜʰˡ)/sum(Δλ)
end

"
    EKoverE = (PCⱼᵐ/(Chl¨C*a̅ⱼ))/((aⱼ*I)/(a̅ⱼ))

Where: 
PCᵐⱼ = maximum carbon-specific growth rate for plankton j
Chl¨Cⱼ = chlorophyll to carbon ratio
aⱼ = slope of irradiance curve
a̅ⱼ = summed slope of irradiance curve (across wavelenghts)
I = irradiance

"
function EKoverE(PCⱼᵐ, Chl¨C, a̅ⱼ, aⱼ, I)
    (PCⱼᵐ/(Chl¨C*a̅ⱼ))/((aⱼ*I)/(a̅ⱼ))
end


"
    if EKoverE<=1
        γⁱⁿʰⱼ = cⁱⁿʰⱼ
    else
        γⁱⁿʰⱼ = 1

Light limitation (Geider formulation)

Where: 
cⁱⁿʰⱼ =
EKoverE = 

"
function γⱼⁱⁿʰ(EKoverE, cⱼⁱⁿʰ) 
    if EKoverE<=1
        return cⱼⁱⁿʰ
    else
        return 1
    end
end

"
    fᵃᴵⱼ = sum(aᶜʰˡⱼₗ, Iₗ)

Summed response to irradiance for all wavelenghts of plankton j.

Where: 
aᶜʰˡⱼₗ = chlorophyll specific slope of irradiance curve at wavelenght l
Iₗ = irradiance at wavelenght l

"
function fⱼᵃᴵ(aⱼₗᶜʰˡ, I)
    sum(aⱼₗᶜʰˡ, Iₗ)
end

"
    aᶜʰˡⱼₗ = Φₘⱼ * apsᶜʰˡⱼₗ

chlorophyll specific slope of irradiance curve for plankton j.

Where: 
Φₘⱼ = maximum quantum yield
apsᶜʰˡⱼₗ = absorption by PS active pigments

"
function aⱼₗᶜʰˡ(Φₘⱼ, apsⱼₗᶜʰˡ)
    Φₘⱼ * apsⱼₗᶜʰˡ
end

export
    darwin_default_light_limitation
    smith_light_limitation
    acclimated_chl_carbon_ratio
    summed_irradiance_curve
    EKoverE
    γⱼⁱⁿʰ
    fⱼᵃᴵ
    aⱼₗᶜʰˡ
end # module