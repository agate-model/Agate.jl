module Models

# -----------------------------------------------------------------------------
# Model modules
# -----------------------------------------------------------------------------

include("NiPiZD/NiPiZD.jl")

export NiPiZD


# The factory types remain available for internal/advanced usage via fully-qualified
# names (e.g. `Agate.Models.NiPiZD.NiPiZDFactory`).
using .NiPiZD: NiPiZDFactory

end # module
