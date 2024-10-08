"""
Modules related to temperature controls on plankton and biogeochemistry

"""
module Temperature

"
    Q10 ^ (T / 10)

Q₁₀ formulation of temperature sensitivity.

Where: 
Q10 = temperature coefficient (value is usually ~2 to 3)
T = Temperature (in degree C or K)

"
function Q₁₀_temperature(Q₁₀, T)
    return Q₁₀^(T / 10)
end

export Q₁₀_temperature
end # module
