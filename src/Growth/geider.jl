"
    PCᵐⱼ = PCᵐᵃˣⱼ * γⁿᵘᵗⱼ * fᵗᵉᵐᵖⱼ *  γᶜᵒ²ⱼ

Realised maximum photosynthetic rate for plankton j.
(ignoring light limitation)

Where: 
PCᵐᵃˣⱼ = maximum photosynthetic rate
γⁿᵘᵗⱼ = nutrient limitation
fᵗᵉᵐᵖⱼ = temperature dependency
γᶜᵒ² = CO2 limitation

"
function PCᵐⱼ(PCᵐᵃˣⱼ, γⁿᵘᵗⱼ, fᵗᵉᵐᵖⱼ, γᶜᵒ²ⱼ)
    PCᵐᵃˣⱼ * γⁿᵘᵗⱼ * fᵗᵉᵐᵖⱼ *  γᶜᵒ²ⱼ
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
function Chl¨Cᵃᶜᶜˡⁱᵐⱼ(Chl¨Cᵐⁱⁿⱼ, Chl¨Cᵐᵃˣⱼ, aⱼ, I, PCᵐⱼ)
    Chl¨Cᵃᶜᶜˡⁱᵐⱼ = (Chl¨Cᵐᵃˣⱼ/((1 + Chl¨Cᵐᵃˣⱼ*aⱼ*I)/(2*PCᵐⱼ)))
    if PCᵐⱼ == 0
        Chl¨Cᵐⁱⁿ = 0
    end
    if Chl¨Cᵃᶜᶜˡⁱᵐ < Chl¨Cᵐⁱⁿⱼ
        return Chl¨Cᵐⁱⁿⱼ
    elseif Chl¨Cᵃᶜᶜˡⁱᵐ > Chl¨Cᵐᵃˣⱼ
        return Chl¨Cᵐᵃˣⱼ
    else
        return Chl¨Cᵃᶜᶜˡⁱᵐⱼ
    end
end

"
    Chl¨Cⱼ = Chl¨Cᵃᶜˡⱼ

Chlorophyll to carbon ratio of plankton j. 
This function can be used to define Chl¨C if variable Chl to Carbon ratios are not used. 

Where: 
Chl¨Cᵃᶜˡ = chlorophyll to carbon ratio based on 
acclimation state of plankton j.


"
function Chl¨C(Chl¨Cᵃᶜˡⱼ)
    Chl¨Cᵃᶜˡⱼ
end

"
    a̅ⱼ = sum(Δλαᶜʰˡⱼₗ)/sum(Δλ)

Summed slope of irradiance curve (across wavelenghts) for plankton j.

Where: 
Δ = 
λ = 
αᶜʰˡⱼₗ = chlorophyll specific slope of irradiance curve

"
function a̅ⱼ()
    sum(Δλαᶜʰˡⱼₗ)/sum(Δλ)
end

"
    EKoverE = (PCᵐⱼ/(Chl¨C*a̅ⱼ))/((aⱼ*I)/(a̅ⱼ))

Where: 
PCᵐⱼ = maximum carbon-specific growth rate for plankton j
Chl¨Cⱼ = chlorophyll to carbon ratio
aⱼ = slope of irradiance curve
a̅ⱼ = summed slope of irradiance curve (across wavelenghts)
I = irradiance

"
function EKoverE(PCᵐⱼ, Chl¨C, a̅ⱼ, aⱼ, I)
    (PCᵐⱼ/(Chl¨C*a̅ⱼ))/((aⱼ*I)/(a̅ⱼ))
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
function γⁱⁿʰⱼ(EKoverE, cⁱⁿʰⱼ) 
    if EKoverE<=1
        return cⁱⁿʰⱼ
    else
        return 1
    end
end

"
    if Iₜₒₜ > Iₘᵢₙ
        PCⱼ = PCᵐⱼ*(1-exp(-(γᶜᶠᵉⱼ*αⱼ*I*Chl¨Cⱼ)/(PCᵐⱼ)))*γⁱⁿʰⱼ
    else
        PCⱼ = 0

Realised carbon-specific growth rate for plankton j (including irradiance).

Where: 
Iₜₒₜ = total irradiance accross wavelenghts
Iₘᵢₙ = minimum irradiance for photosynthesis to occur 
PCᵐⱼ = realised maximum growth rate (ignoring light) of plankton j
γᶜᶠᵉⱼ = effect of iron limitation on photosynthesis of plankton j
Chl¨Cⱼ = chlorophyll to carbon ratio of plankton j
γⁱⁿʰⱼ = light inhibition of plankton j
fᵃᴵⱼ = summed response to irradiance for all wavelenghts of plankton j

"
function geider_PCⱼ(PCᵐⱼ,  γⁱⁿʰⱼ, Iₜₒₜ, Iₘᵢₙ, Chl¨Cⱼ)
    if Iₜₒₜ > Iₘᵢₙ
        return PCᵐⱼ*(1-exp(-(γᶜᶠᵉ*fᵃᴵⱼ*Chl¨Cⱼ)/(PCᵐⱼ)))*γⁱⁿʰⱼ
    else
        return 0
    end
end

"
    fᵃᴵⱼ = sum(aᶜʰˡⱼₗ, Iₗ)

Summed response to irradiance for all wavelenghts of plankton j.

Where: 
aᶜʰˡⱼₗ = chlorophyll specific slope of irradiance curve at wavelenght l
Iₗ = irradiance at wavelenght l

"
function fᵃᴵⱼ(aᶜʰˡⱼₗ, I)
    sum(aᶜʰˡⱼₗ, Iₗ)
end

"
    aᶜʰˡⱼₗ = Φₘⱼ * apsᶜʰˡⱼₗ

chlorophyll specific slope of irradiance curve for plankton j.

Where: 
Φₘⱼ = maximum quantum yield
apsᶜʰˡⱼₗ = absorption by PS active pigments

"
function aᶜʰˡⱼₗ(Φₘⱼ, apsᶜʰˡⱼₗ)
    Φₘⱼ * apsᶜʰˡⱼₗ
end