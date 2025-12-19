module NiPiZD

include("Tracers.jl")
include("Constructor.jl")

using .Tracers
using .Constructor

export construct
export instantiate

export default_phyto_pft_parameters
export default_phyto_geider_pft_parameters
export default_zoo_pft_parameters
export default_bgc_specification

end # module
