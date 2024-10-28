"""
Modules related to photosynthetically available radiation (PAR)
"""

module Light

using Oceananigans
using Oceananigans: Clock
using Oceananigans.Units
using Oceananigans.Fields: FunctionField, compute!

import Oceananigans.Biogeochemistry:
    update_biogeochemical_state!, biogeochemical_auxiliary_fields

const year = years = 365day

"""
    cyclical_PAR(t, z) -> Float

Time-dependent cyclical PAR at depth `z` (suitable for use with box models).
"""
function cyclical_PAR(z, t)
    PAR⁰ =
        60 *
        (1 - cos((t + 15days) * 2π / year)) *
        (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
    return PAR⁰ * exp(0.2 * z)
end

cyclical_PAR(; z) = t -> cyclical_PAR(z, t)

"""
Time-dependent cyclical PAR at depth `z` (suitable for use with column models).
"""
function cyclical_PAR(x, y, z, t)
    return cyclical_PAR(z, t)
end

"""
Light module for cyclical PAR (can be used with box or column models).

# Fields
- `field`: Oceananigans.FunctionField
"""
struct CyclicalPhotosyntheticallyActiveRadiation
    field
end

"""
    CyclicalPhotosyntheticallyActiveRadiation() -> DataType

# Keywords
- `grid`: the geometry to build the model on defined as an Oceananigans grid object
- `PAR_f`: a time dependant PAR function (defaults to `Agate.Library.Light.cyclical_PAR`)
"""
function CyclicalPhotosyntheticallyActiveRadiation(; grid, PAR_f=cyclical_PAR(; z=-10))
    clock = Clock(; time=zero(grid))
    PAR = FunctionField{Center,Center,Center}(PAR_f, grid; clock)
    return CyclicalPhotosyntheticallyActiveRadiation(PAR)
end

function update_biogeochemical_state!(model, PAR::CyclicalPhotosyntheticallyActiveRadiation)
    PAR.field.clock.time = model.clock.time
    compute!(PAR.field)
    return nothing
end

function biogeochemical_auxiliary_fields(par::CyclicalPhotosyntheticallyActiveRadiation)
    return (PAR=par.field,)
end

export cyclical_PAR, CyclicalPhotosyntheticallyActiveRadiation

end # module
