module Models

include("Parameters.jl")
include("NiPiZD/NiPiZD.jl")
include("DARWIN/DARWIN.jl")

using .Parameters
using .NiPiZD
using .DARWIN

export Parameters
export NiPiZD
export DARWIN

export AbstractDiameterSpecification
export DiameterListSpecification
export DiameterRangeSpecification
export NiPiZDBiogeochemistrySpecification
export DarwinBiogeochemistrySpecification
export PhytoPFTParameters
export ZooPFTParameters
export PhytoSpecification
export ZooSpecification
export NiPiZDParameters
export DarwinParameters
export create_nipizd_parameters
export compute_nipizd_parameters
export create_darwin_parameters
export compute_darwin_parameters

end # module
