module Library

# #Growth
# export
#     default_PCⱼ

# #Light    
# export 
#     γⱼˡⁱᵍʰᵗ
#     smith_light_limitation

# # Nutrients
# export 
#     monod_limitation

# #Temperature
# export
#     Q₁₀_temperature


include("Allometry/Allometry.jl")
include("Chlorophyll/Chlorophyll.jl")
include("Growth/Growth.jl")
include("Light/Light.jl")
include("Mortality/Mortality.jl")
include("Nutrients/Nutrients.jl")
include("Predation/Predation.jl")
include("Remineralization/Remineralization.jl")
include("Temperature/Temperature.jl")

using .Allometry
using .Chlorophyll
using .Growth
using .Light
using .Mortality
using .Nutrients
using .Predation
using .Remineralization
using .Temperature

end # module