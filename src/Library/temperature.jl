"""Temperature-response kernels."""

module Temperature

export q10_temperature_factor

"""
    q10_temperature_factor(T, Q10)
    q10_temperature_factor(T, Q10, reference_temperature)

Compute the Q10 temperature factor relative to a reference temperature. The
two-argument form uses a reference temperature of zero.

# Arguments
- `T`: temperature in degrees Celsius
- `Q10`: Q10 coefficient
- `reference_temperature`: reference temperature in degrees Celsius
"""
@inline q10_temperature_factor(T, Q10, reference_temperature) =
    Q10^((T - reference_temperature) / 10)
@inline q10_temperature_factor(T, Q10) = q10_temperature_factor(T, Q10, zero(T))

end # module
