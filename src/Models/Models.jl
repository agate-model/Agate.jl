"""Bundled model-family implementations."""
module Models

# -----------------------------------------------------------------------------
# Model modules
# -----------------------------------------------------------------------------

include("NiPiZD/NiPiZD.jl")
include("FrankenLOBSTER/FrankenLOBSTER.jl")

export NiPiZD
export FrankenLOBSTER

end # module
