"""
A library of modules to create marine biogeochemical models

"""
module Library

include("allometry.jl")
include("chlorophyll.jl")
include("growth.jl")
include("photosynthesis.jl")
include("mortality.jl")
include("nutrients.jl")
include("predation.jl")
include("remineralization.jl")
include("temperature.jl")

using .Allometry
using .Chlorophyll
using .Growth
using .Mortality
using .Nutrients
using .Photosynthesis
using .Predation
using .Remineralization
using .Temperature

end # module
