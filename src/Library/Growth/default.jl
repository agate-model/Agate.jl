"
    PCⱼ = PCⱼᵐᵃˣ * γⱼⁿᵘᵗ *  γⱼˡⁱᵍʰᵗ * fⱼᵗᵉᵐᵖ *  γⱼᶜᵒ²

Carbon-specific growth rate for plankton j (Default MITgcm-DARWIN formulation).

Where: 
PCᵐᵃˣⱼ = maximum carbon-specific growth rate for plankton j,
γⁿᵘᵗⱼ = nutrient limition,
γˡⁱᵍʰᵗⱼ = light limition,
fᵗᵉᵐᵖⱼ = temperature limitation,
γᶜᵒ²ⱼ = carbon dioxide limitation

"
function default_PCⱼ(PCⱼᵐᵃˣ, γⱼⁿᵘᵗ,  γⱼˡⁱᵍʰᵗ, fⱼᵗᵉᵐᵖ,  γⱼᶜᵒ²)
    PCⱼᵐᵃˣ * γⱼⁿᵘᵗ *  γⱼˡⁱᵍʰᵗ * fⱼᵗᵉᵐᵖ *  γⱼᶜᵒ²
end