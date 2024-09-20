module Library

export
    Allometry
    Chlorophyll
    Growth
    Light
    Mortality
    Nutrients
    Predation
    Remineralization
    Temperature

include("Allometry/Allometry.jl")
include("Chlorophyll/Chlorophyll.jl")
include("Growth/Growth.jl")
include("Light/Light.jl")
include("Mortality/Mortality.jl")
include("Nutrients/Nutrients.jl")
include("Predation/Predation.jl")
include("Remineralization/Remineralization.jl")
include("Temperature/Temperature.jl")

end # module