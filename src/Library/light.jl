"""Photosynthetically Active Radiation (PAR) utilities."""

module Light

using Oceananigans.Units

const year = years = 365day

export CyclicalPAR, cyclical_par_at_depth

"""
    cyclical_par_at_depth(z, t)

Evaluate the idealized seasonal photosynthetically active radiation (PAR) field
at depth `z` and time `t`.

!!! formulation
    ```math
    I_0(t) = 60\\left[1 - \\cos\\left(\\frac{2\\pi(t + 15\\,days)}{year}\\right)\\right]
             \\left[1 + 0.2\\exp\\left(-\\left(\\frac{mod(t, year)-200\\,days}{50\\,days}\\right)^2\\right)\\right]^{-1} + 2
    ```

    with vertical attenuation

    ```math
    I(z,t) = I_0(t)\\exp(0.2z).
    ```

    In Oceananigans coordinates, negative `z` is below the surface, so
    `exp(0.2z)` attenuates light with depth.
"""
@inline function cyclical_par_at_depth(z, t)
    PAR⁰ =
        60 *
        (1 - cos((t + 15days) * 2π / year)) *
        (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
    return PAR⁰ * exp(0.2 * z)
end

"""
    CyclicalPAR(z)

Cyclical, depth-attenuated PAR evaluated at fixed depth `z`.
`CyclicalPAR(z)(t)` is convenient for box models and direct time-series evaluation.
For gridded models, use `cyclical_par_at_depth(z, t)` in a spatial function passed to
`Oceananigans.Fields.FunctionField`.

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

@inline (p::CyclicalPAR)(t) = cyclical_par_at_depth(p.z, t)

end # module
