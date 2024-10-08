"""
Modules related to photosynthetically available radiation (PAR)

"""

module Light

using Oceananigans.Units

const year = years = 365day

"""
    cyclical_PAR(t, z) -> Float

Time-dependent cyclical PAR at depth z (suitable for use with box models).
"""
function cyclical_PAR(t; z=-10)
    PAR⁰ =
        60 *
        (1 - cos((t + 15days) * 2π / year)) *
        (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
    return PAR⁰ * exp(0.2z)
end

export cyclical_PAR

end # module
