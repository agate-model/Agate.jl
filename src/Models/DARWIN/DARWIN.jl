module DARWIN

include("Tracers.jl")
include("DarwinParameters.jl")
include("Constructor.jl")

using .DarwinParameters
using .Tracers
using .Constructor

export construct
export instantiate

export default_phyto_pft_parameters
export default_phyto_geider_pft_parameters
export default_zoo_pft_parameters
export default_darwin_bgc_specification

end # module
