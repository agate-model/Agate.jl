"""Light-response kernels used by phytoplankton growth formulations."""
module Photosynthesis

export smith_light_limitation, geider_light_response, exponential_light_limitation

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
    exponential_light_limitation(PAR, light_scale)

Evaluate a saturating exponential light-response factor.

```math
L_I(I) = 1 - \\exp\\left(-\\frac{I}{K_I}\\right)
```

`PAR` is photosynthetically active radiation and `light_scale` is the positive
irradiance scale ``K_I``. A documented ocean-biogeochemical use of this response is
Lévy, Klein & Tréguier (2001), Eq. (A7), *Journal of Marine Research* 59, 535-565,
doi:10.1357/002224001762842181.
"""
@inline exponential_light_limitation(PAR, light_scale) =
    one(PAR + light_scale) - exp(-PAR / light_scale)


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
