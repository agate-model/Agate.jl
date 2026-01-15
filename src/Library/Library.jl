"""
A library of modules to create marine biogeochemical models

"""
module Library

# Construction-time symbolic equation system (no Models/Constructor deps).
include("Equations.jl")

include("allometry.jl")
include("light.jl")
include("mortality.jl")
include("nutrients.jl")
include("photosynthesis.jl")
include("predation.jl")
include("remineralization.jl")
include("temperature.jl")

using .Allometry
using .Equations
using .Light
using .Mortality
using .Nutrients
using .Photosynthesis
using .Predation
using .Remineralization
using .Temperature

end # module
