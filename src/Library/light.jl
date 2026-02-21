"""Photosynthetically Active Radiation (PAR) utilities."""

module Light

using Adapt
using Oceananigans
using Oceananigans: Clock
using Oceananigans.Units
using Oceananigans.Fields: FunctionField, compute!

import Oceananigans.Biogeochemistry:
    update_biogeochemical_state!, biogeochemical_auxiliary_fields

const year = years = 365day

export CyclicalPAR, FunctionFieldPAR

@inline function _cyclical_par_at_depth(z, t)
    PAR⁰ =
        60 *
        (1 - cos((t + 15days) * 2π / year)) *
        (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
    return PAR⁰ * exp(0.2 * z)
end

"""
    CyclicalPAR()
    CyclicalPAR(z)

Cyclical, depth-attenuated PAR suitable for box and column models.

- For box models: `CyclicalPAR(z)(t)` evaluates PAR at fixed depth `z`.
- For column models: `CyclicalPAR()(x, y, z, t)` evaluates PAR at runtime depth `z`.

!!! formulation
    60 *
    (1 - cos((t + 15days) * 2π / year)) *
    (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

    with depth attenuation:

    PAR(z, t) = PAR⁰(t) * exp(0.2 * z)
"""
struct CyclicalPAR{Z}
    z::Z
end

@inline CyclicalPAR() = CyclicalPAR(nothing)

@inline (p::CyclicalPAR{T})(t) where {T} = _cyclical_par_at_depth(p.z, t)

@inline function (p::CyclicalPAR{Nothing})(x, y, z, t)
    return _cyclical_par_at_depth(z, t)
end

"""
    FunctionFieldPAR(field)

Light module wrapping an `Oceananigans.FunctionField` that represents PAR.

# Fields
- `field`: `Oceananigans.Fields.FunctionField`
"""
struct FunctionFieldPAR{F}
    field::F
end

Adapt.@adapt_structure FunctionFieldPAR

"""
    FunctionFieldPAR(; grid, PAR_f=CyclicalPAR(-10))

Create a PAR module backed by a `FunctionField`.

# Keywords
- `grid`: an Oceananigans grid.
- `PAR_f`: a callable compatible with `FunctionField` (defaults to `CyclicalPAR(-10)`).
"""
function FunctionFieldPAR(; grid, PAR_f=CyclicalPAR(-10))
    clock = Clock(; time=0.0)
    PAR_field = FunctionField{Center,Center,Center}(PAR_f, grid; clock)
    return FunctionFieldPAR(PAR_field)
end

"""
    update_biogeochemical_state!(model, PAR::FunctionFieldPAR)

`Oceananigans.Biogechemistry` hook that computes and updates the irradiance field in-place.
"""
function update_biogeochemical_state!(model, PAR::FunctionFieldPAR)
    PAR.field.clock.time = model.clock.time
    compute!(PAR.field)
    return nothing
end

"""
    biogeochemical_auxiliary_fields(par::FunctionFieldPAR)

Return a named tuple containing the Photosynthetically Active Radiation (PAR) field.
"""
function biogeochemical_auxiliary_fields(par::FunctionFieldPAR)
    return (PAR=par.field,)
end

end # module
