"""
Functions related to photosynthetically available radiation (PAR).
"""

module Light

using Oceananigans
using Oceananigans: Clock
using Oceananigans.Units
using Oceananigans.Fields: FunctionField, compute!

import Oceananigans.Biogeochemistry:
    update_biogeochemical_state!, biogeochemical_auxiliary_fields

const year = years = 365day

export cyclical_PAR, FunctionFieldPAR

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
Time and depth dependent cyclical PAR (suitable for use with column models).
"""
function cyclical_PAR(x, y, z, t)
    return cyclical_PAR(z, t)
end

"""
Light module for PAR defined by simple functions (can be used with box or column models).

# Fields
- `field`: Oceananigans.FunctionField
"""
struct FunctionFieldPAR
    field
end

"""
    FunctionFieldPAR() -> DataType

# Keywords
- `grid`: the geometry to build the model on defined as an Oceananigans grid object
- `PAR_f`: a PAR function of time (and depth), defaults to `cyclical_PAR`
"""
function FunctionFieldPAR(; grid, PAR_f=cyclical_PAR(; z=-10))
    clock = Clock(; time=zero(grid))
    PAR_field = FunctionField{Center,Center,Center}(PAR_f, grid; clock)
    return FunctionFieldPAR(PAR_field)
end

function update_biogeochemical_state!(model, PAR::FunctionFieldPAR)
    PAR.field.clock.time = model.clock.time
    compute!(PAR.field)
    return nothing
end

function biogeochemical_auxiliary_fields(par::FunctionFieldPAR)
    return (PAR=par.field,)
end

end # module
