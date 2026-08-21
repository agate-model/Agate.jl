"""Temperature-response formulations."""

module Temperature

export q10_temperature_factor

"""
    q10_temperature_factor(T, Q10)

Compute the Q10 temperature factor.

!!! formulation
    Q10 ^ (T / 10)

# Arguments
- `T`: temperature (°C)
- `Q10`: Q10 coefficient
"""
@inline q10_temperature_factor(T, Q10) = Q10^(T / 10)

end # module
