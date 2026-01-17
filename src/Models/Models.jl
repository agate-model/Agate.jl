module Models

using ..FactoryInterface:
    factory_groups,
    default_plankton_dynamics,
    default_community,
    default_biogeochem_dynamics

include("InteractionDefaults.jl")
using .InteractionDefaults: default_palatability_provider, default_assimilation_provider
export default_palatability_provider, default_assimilation_provider


# -----------------------------------------------------------------------------
# Model modules
# -----------------------------------------------------------------------------

include("NiPiZD/NiPiZD.jl")
include("DARWIN/DARWIN.jl")

export NiPiZD, DARWIN

# The factory types remain available for internal/advanced usage via fully-qualified
# names (e.g. `Agate.Models.NiPiZD.NiPiZDFactory`).
using .NiPiZD: NiPiZDFactory
using .DARWIN: DarwinFactory

end # module
