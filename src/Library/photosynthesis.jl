"""Light-response kernels used by phytoplankton growth formulations."""
module Photosynthesis

export smith_light_limitation, geider_light_response

"""
    smith_light_limitation(PAR, alpha, maximum_rate)

Evaluate the dimensionless Smith (1936) light-limitation factor.

```math
L_S(I) = \\frac{\\alpha I}{\\sqrt{\\mu_{max}^2 + (\\alpha I)^2}}
```

`PAR` is photosynthetically active radiation, `alpha` is the initial photosynthetic
slope, and `maximum_rate` is the enclosing growth-process rate scale.
"""
@inline function smith_light_limitation(PAR, alpha, maximum_rate)
    if alpha == zero(alpha) || maximum_rate == zero(maximum_rate)
        return zero(alpha)
    end
    light_rate = alpha * PAR
    return light_rate / sqrt(maximum_rate * maximum_rate + light_rate * light_rate)
end

"""
    geider_light_response(PAR, alpha, maximum_rate, chlorophyll_to_carbon_ratio)

Evaluate the dimensionless Geider light-response factor.

```math
L_G(I) = 1 - \\exp\\left(-\\frac{\\alpha^{chl}\\theta^C I}{\\mu_{max}}\\right)
```

`PAR` is photosynthetically active radiation, `alpha` is the chlorophyll-specific
initial slope, `chlorophyll_to_carbon_ratio` is ``\\theta^C``, and `maximum_rate` is
the enclosing growth-process rate scale. Multiplying the returned factor by
`maximum_rate` gives the light-dependent growth scale.
"""
@inline function geider_light_response(
    PAR, alpha, maximum_rate, chlorophyll_to_carbon_ratio
)
    maximum_rate == zero(maximum_rate) && return zero(maximum_rate)
    return one(maximum_rate) - exp(
        (-alpha * chlorophyll_to_carbon_ratio * PAR) / maximum_rate
    )
end

end # module
