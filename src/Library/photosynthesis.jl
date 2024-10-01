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
μ₀ = Maximum growth rate at T = 0 °C

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
    fᵃᴵⱼ = aⱼᶜʰˡ * Iₜₒₜ

Non-spectral maximum carbon yield of photosynthesis.

In this formulation, the maximum carbon yield is estimated by using
a average spectral absorption co-efficient and then multiplying this
by the total PAR. 

Where: 
aᶜʰˡⱼ = average chlorophyll specific slope of irradiance curve for plankton j,
Iₜₒₜ = total PAR.

"
function non_spectral_carbon_yield(aⱼᶜʰˡ, Iₜₒₜ)
    aⱼᶜʰˡ * Iₜₒₜ
end

"
    fᵃᴵⱼ = sum(aᶜʰˡⱼₗ, Iₗ)

Spectral maximum carbon yield of photosynthesis based on all wavelenghts for 
plankton j.

In this function, the total carbon yield is estimated by summing the 
carbon yield for each spectrum based on the coefficient of absorption 
by photosynthetically active pigments.

Where: 
aᶜʰˡⱼₗ = chlorophyll specific slope of irradiance curve at wavelenght l
Iₗ = irradiance at wavelenght l

"
function spectral_carbon_yield(aⱼₗᶜʰˡ, Iₗ)
    sum(aⱼₗᶜʰˡ, Iₗ)
end

"
    Chl¨Cⱼ = Chlⱼ/Cⱼ

Estimates chlorophyll to carbon ratio based on tracer values.

"
function chlorophyll_carbon_ratio(Chlⱼ, Cⱼ)
    Chlⱼ/Cⱼ
end

"
    γⱼˡⁱᵐ = (1-exp(-(γⱼᶜᶠᵉ*fⱼᵃᴵ*Chl¨Cⱼ)/(PCᵐⱼ)))

Estimates geider light limitation 

Where: 
Chl¨Cⱼ = chlorophyll to carbon ratio
γⱼᶜᶠᵉ = iron limitation
fⱼᵃᴵ = carbon yield from photosynthesis
PCᵐⱼ = geider light saturated growth

"
function geider_light_limitation(Chl¨Cⱼ, γⱼᶜᶠᵉ, fⱼᵃᴵ, PCᵐⱼ)
    (1-exp(-(γⱼᶜᶠᵉ*fⱼᵃᴵ*Chl¨Cⱼ)/(PCᵐⱼ)))
end

"
    Chl¨Cᵃᶜˡⱼ = (Chl¨Cᵐᵃˣⱼ/((1 + Chl¨Cᵐᵃˣⱼ*aⱼ*I)/(2*PCᵐⱼ)))

Acclimated chlorophyll quota for plankton j 
(computed indepently of chlorophyll and carbon tracers).

Where: 
Chl¨Cᵐⁱⁿⱼ = min Chl to Carbon ratio,
Chl¨Cᵐᵃˣⱼ = max Chl to Carbon ratio,
I = irradiance,
aⱼ = slope of irradiance curve,
PCᵐⱼ = realised maximum growth rate (ignoring light).

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
    a̅ⱼ = sum(Δλₗ*αⱼₗᶜʰˡ)/sum(Δλ)

Summed slope of irradiance curve (across wavelenghts) for plankton j.

Where: 
Δ = delta wavelenght 
λ = wavelenght
αᶜʰˡⱼₗ = chlorophyll specific slope of irradiance curve at wavelenght l

"
function summed_irradiance_curve(Δλₗ, αⱼₗᶜʰˡ)
    sum(Δλₗ*αⱼₗᶜʰˡ)/sum(Δλₗ)
end

"
    aᶜʰˡⱼₗ = Φₘⱼ * apsᶜʰˡⱼₗ

chlorophyll-a specific slope of irradiance curve for plankton j.

Where: 
Φₘⱼ = maximum quantum yield
apsᶜʰˡⱼₗ = absorption by PS active pigments

"
function chlorophyll_irradiance_slope(Φₘⱼ, apsⱼₗᶜʰˡ)
    Φₘⱼ * apsⱼₗᶜʰˡ
end

"
    if EKoverE<=1
        γⁱⁿʰⱼ = cⁱⁿʰⱼ
    else
        γⁱⁿʰⱼ = 1

Photo inhibition (Geider formulation)

Where: 
cⁱⁿʰⱼ = photo-inhibition coefficient for Geider growth
EKoverE = growth potential based on chl:c and light response

"
function geider_light_inhibition(EKoverE, cⱼⁱⁿʰ) 
    if EKoverE<=1
        return cⱼⁱⁿʰ
    else
        return 1
    end
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

export
    smith_light_limitation
    darwin_default_light_limitation
    geider_light_limitation
    non_spectral_carbon_yield
    spectral_carbon_yield
    chlorophyll_carbon_ratio
    acclimated_chl_carbon_ratio
    summed_irradiance_curve
    chlorophyll_irradiance_slope
    geider_light_inhibition
    EKoverE
end # module