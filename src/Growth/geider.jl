"
    if Iₜₒₜ > Iₘᵢₙ
        (1-exp(-(γᶜᶠᵉⱼ*αⱼ*I*Chl¨Cⱼ)/(PCᵐⱼ)))*γⁱⁿʰⱼ
    else
        PCⱼ = 0

Light limitation for plankton j (Geider MITgcm-DARWIN formulation). 

Where: 
Iₜₒₜ = total irradiance accross wavelenghts
Iₘᵢₙ = minimum irradiance for photosynthesis to occur 
γᶜᶠᵉⱼ = effect of iron limitation on photosynthesis of plankton j
Chl¨Cⱼ = chlorophyll to carbon ratio of plankton j
γⁱⁿʰⱼ = light inhibition of plankton j
fᵃᴵⱼ = summed response to irradiance for all wavelenghts of plankton j

"
function geider_γⱼˡⁱᵍʰᵗ(γⱼⁱⁿʰ, Iₜₒₜ, Iₘᵢₙ, Chl¨Cⱼ, γⱼᶜᶠᵉ, fⱼᵃᴵ)
    if Iₜₒₜ > Iₘᵢₙ
        return (1-exp(-(γⱼᶜᶠᵉ*fⱼᵃᴵ*Chl¨Cⱼ)/(PCᵐⱼ)))*γⱼⁱⁿʰ
    else
        return 0
    end
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
function Chl¨Cⱼᵃᶜˡ(Chl¨Cⱼᵐⁱⁿ, Chl¨Cⱼᵐᵃˣ, aⱼ, I, PCⱼᵐ)
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
    Chl¨Cⱼ = Chl¨Cⱼᵃᶜˡ

Chlorophyll to carbon ratio of plankton j. 
This function can be used to define Chl¨C if variable Chl to Carbon ratios are not used. 

Where: 
Chl¨Cⱼᵃᶜˡ = chlorophyll to carbon ratio based on 
acclimation state of plankton j.


"
function Chl¨Cⱼ(Chl¨Cⱼᵃᶜˡ)
    Chl¨Cⱼᵃᶜˡ
end

"
    a̅ⱼ = sum(Δλαⱼₗᶜʰˡ)/sum(Δλ)

Summed slope of irradiance curve (across wavelenghts) for plankton j.

Where: 
Δ = 
λ = 
αᶜʰˡⱼₗ = chlorophyll specific slope of irradiance curve

"
function a̅ⱼ()
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