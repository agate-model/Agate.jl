"""Temperature-response functors."""

module Temperature

export q10_temperature_factor

"""
    Q10Temperature(Q10)

Q₁₀ formulation of temperature sensitivity.

!!! formulation
    Q10 ^ (T / 10)

where:
- `Q10` = temperature coefficient (typically ≈ 2–3)
- `T` = temperature (°C)
"""
struct Q10Temperature{T}
    Q10::T
end

@inline (q::Q10Temperature)(T) = q.Q10^(T / 10)

"""
    q10_temperature_factor(T, Q10)

Compute the Q₁₀ temperature factor.

!!! formulation
    Q10 ^ (T / 10)

This is a thin, inlined alias around `Q10Temperature(Q10)(T)`.

# Arguments
- `T`: temperature (°C)
- `Q10`: Q₁₀ coefficient
"""
@inline q10_temperature_factor(T, Q10) = Q10Temperature(Q10)(T)

end # module
