module Models

include("Parameters.jl")
include("NiPiZD/NiPiZD.jl")

using .Parameters
using .NiPiZD

export Parameters
export NiPiZD

export AbstractDiameterSpecification
export DiameterListSpecification
export DiameterRangeSpecification
export NiPiZDBiogeochemistrySpecification
export PhytoPFTParameters
export ZooPFTParameters
export PhytoSpecification
export ZooSpecification
export NiPiZDParameters
export create_nipizd_parameters
export compute_nipizd_parameters

end # module
