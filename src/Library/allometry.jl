module Allometry

"""
    allometric_scaling_power(a, b, V::Number)

Allometric scaling function using the power law for cell volume.

Note that for cells with a diameter < ~6 micrometers, `b` tends to be negative.
For cells > ~6 micrometers, `b` tends to be positive.

# Arguments
- `a`: scale
- `b`: exponent
- `V`: cell volume
"""
function allometric_scaling_power(a::Number, b::Number, V::Number)
    return a * V^b
end

"""
    allometric_scaling_power(a, b, d::Number)

Allometric scaling function using the power law for cell diameter. Converts
diameter to volume assuming a spherical shape.

Note that if diameter is passed instead of volume this should be done explicitly:
`x = (2.19, -0.16, d=10)`

# Arguments
- `a`: scale
- `b`: exponent
- `d`: cell diameter
"""
function allometric_scaling_power(a::Number, b::Number, d::Number)
    V = (4 / 3) * π * (d / 2)^3
    return allometric_scaling_power(a, b, V)
end

end # module
