"""
Modules related to photosynthetically available radiation (PAR)

"""

module Light

using Oceananigans.Units

const year = years = 365day


"""
    PAR⁰(t) -> Float

Time-dependent surface-level PAR.
"""
function PAR⁰(t)
    return 60 *
    (1 - cos((t + 15days) * 2π / year)) *
    (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
end


"""
    PAR_box(t, z) -> Float

Time-dependent PAR at depth z.
"""
function PAR_box(t; z=-10)
    return PAR⁰(t) * exp(0.2z)
end

export
    PAR_box

end # module
