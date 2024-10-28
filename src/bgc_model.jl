include("Library/light.jl")

using .Light

using OceanBioME
using Oceananigans
using Oceananigans: Clock

using Oceananigans.Units
using Oceananigans.Fields: FunctionField, compute!

using OceanBioME: Biogeochemistry

import Oceananigans: set!
import Oceananigans.Biogeochemistry:
    update_biogeochemical_state!, biogeochemical_auxiliary_fields

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
    bgc_tracers;
    grid=BoxModelGrid(),
    light_attenuation=BoxPhotosyntheticallyActiveRadiation(; grid=grid),
)
    return Biogeochemistry(bgc_tracers; light_attenuation=light_attenuation)
end

"""
Light module for use with Box models

# Fields
- `field`: Oceananigans.FunctionField
"""
struct BoxPhotosyntheticallyActiveRadiation
    field
end

"""
    BoxPhotosyntheticallyActiveRadiation() -> DataType

# Keywords
- `grid`: the geometry to build the model on defined as an Oceananigans grid object
- `PAR_f`: a time dependant PAR function (defaults to `Agate.Library.Light.cyclical_PAR`)
"""
function BoxPhotosyntheticallyActiveRadiation(;
    grid=BoxModelGrid(), PAR_f=cyclical_PAR(; z=-10)
)
    clock = Clock(; time=zero(grid))
    PAR = FunctionField{Center,Center,Center}(PAR_f, grid; clock)
    return BoxPhotosyntheticallyActiveRadiation(PAR)
end

function update_biogeochemical_state!(model, PAR::BoxPhotosyntheticallyActiveRadiation)
    PAR.field.clock.time = model.clock.time
    compute!(PAR.field)
    return nothing
end

function biogeochemical_auxiliary_fields(par::BoxPhotosyntheticallyActiveRadiation)
    return (PAR=par.field,)
end
