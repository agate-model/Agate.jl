include("Library/light.jl")

using .Light

using OceanBioME
using OceanBioME: Biogeochemistry

export create_bgc_model

"""
# Arguments
- `bgc_tracers`: biogeochemistry model tracers, a subtype of AbstractContinuousFormBiogeochemistry
    e.g., returned by `Agate.define_tracer_functions`

# Keywords
- `grid`: the geometry to build the model on defined as an Oceananigans grid object
- `light_attenuation`: model for attenuation of PAR through water
"""
function create_bgc_model(
    bgc_tracers; grid=BoxModelGrid(), light_attenuation=FunctionPAR(; grid=grid)
)
    return Biogeochemistry(bgc_tracers; light_attenuation=light_attenuation)
end
