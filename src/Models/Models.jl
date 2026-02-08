module Models

using ..Factories: default_plankton_dynamics, default_community, default_biogeochem_dynamics

# -----------------------------------------------------------------------------
# Shared model utilities
# -----------------------------------------------------------------------------

include("Sums.jl")
include("Variants.jl")
include("InteractionMatrixDerivations.jl")

# -----------------------------------------------------------------------------
# Model modules
# -----------------------------------------------------------------------------

include("NiPiZD/NiPiZD.jl")
include("DARWIN/DARWIN.jl")

export NiPiZD, DARWIN

# Variant scaffolding (developer-oriented).
export ModelId, VariantSpec, register_variant, variant, list_variants

# The factory types remain available for internal/advanced usage via fully-qualified
# names (e.g. `Agate.Models.NiPiZD.NiPiZDFactory`).
using .NiPiZD: NiPiZDFactory
using .DARWIN: DarwinFactory

end # module
