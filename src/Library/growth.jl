"""
Modules related to phytoplankton photosynthetic growth

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
γᶜᵒ²ⱼ = carbon dioxide limitation

"
function default_PCⱼ(PCⱼᵐᵃˣ, γⱼⁿᵘᵗ, γⱼˡⁱᵍʰᵗ, fⱼᵗᵉᵐᵖ, γⱼᶜᵒ²)
    return PCⱼᵐᵃˣ * γⱼⁿᵘᵗ * γⱼˡⁱᵍʰᵗ * fⱼᵗᵉᵐᵖ * γⱼᶜᵒ²
end

export default_PCⱼ
end # module
