"""
Functions related to phytoplankton photosynthetic growth.
"""

module Growth

"
    PC = PCᵐᵃˣ * γⁿᵘᵗ *  γˡⁱᵍʰᵗ * fᵗᵉᵐᵖ *  γᶜᵒ²

Carbon-specific growth rate for plankton (Default MITgcm-DARWIN formulation).

# Arguments
- `PCᵐᵃˣ`: maximum carbon-specific growth rate
- `γⁿᵘᵗ`: nutrient limition
- `γˡⁱᵍʰᵗ`: light limition
- `fᵗᵉᵐᵖ`: temperature limitation
- `γᶜᵒ²`: carbon dioxide limitation
"
function default_PC(PCᵐᵃˣ, γⁿᵘᵗ, γˡⁱᵍʰᵗ, fᵗᵉᵐᵖ, γᶜᵒ²)
    return PCᵐᵃˣ * γⁿᵘᵗ * γˡⁱᵍʰᵗ * fᵗᵉᵐᵖ * γᶜᵒ²
end

export default_PC
end # module
