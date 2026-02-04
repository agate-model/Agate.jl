"""Temperature-response functors."""

module Temperature

export q10_temperature_factor

"""
    Q10Temperature(Q10)

Temperature sensitivity factor `T ↦ Q10^(T/10)`.

`Q10` is the multiplicative rate increase for a 10° temperature rise.
"""
struct Q10Temperature{T}
    Q10::T
end

@inline (q::Q10Temperature)(T) = q.Q10^(T / 10)

"""Q10 temperature factor `Q10^(T/10)`.

This is a thin, inlined alias around `Q10Temperature(Q10)(T)`.
"""
@inline q10_temperature_factor(T, Q10) = Q10Temperature(Q10)(T)

end # module
