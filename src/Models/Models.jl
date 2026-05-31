module Models

# -----------------------------------------------------------------------------
# Shared model utilities
# -----------------------------------------------------------------------------

include("sums.jl")
include("manifests.jl")
# -----------------------------------------------------------------------------
# Model modules
# -----------------------------------------------------------------------------

include("NiPiZD/NiPiZD.jl")
include("DARWIN/DARWIN.jl")

export NiPiZD, DARWIN

export export_manifest, construct_from_manifest

# The factory types remain available for internal/advanced usage via fully-qualified
# names (e.g. `Agate.Models.NiPiZD.NiPiZDFactory`).
using .NiPiZD: NiPiZDFactory
using .DARWIN: DarwinFactory

end # module
