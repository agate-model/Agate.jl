module Mortality

export linear_loss,
    quadratic_loss,
    net_linear_loss_quota,
    net_quadratic_loss_quota

"""
Linear mortality rate.
In this formulation mortality is constant, and can be interpreted as
a "closure term" for low density predation and and other death terms.

# Arguments
- `P`: plankton concentration
- `l`: mortality rate
"""
linear_loss(P, l) = l * P

"""
Quadratic mortality coefficient.
In this formulation mortality increases exponentially with plankton biomass
and is often interpreted to represent viral processes and non-represented density-dependent predation effects.

# Arguments
- `P`: plankton concentration
- `l`: mortality rate
"""
quadratic_loss(P, l) = l * P^2

"""
    net_linear_loss_quota(P, linear_mortality, DOM_POM_fractionation, quota)

Net loss of all plankton due to linear mortality with a elemental quota term.
The quota term is multiplied with estimated loss to convert from e.g. carbon to nitrogen.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `linear_mortality`: NamedArray of all plankton linear mortality rates
- `DOM_POM_fractionation`: Float which represents dissolved-to-particulate fractionation
- `quota`: Float which represents plankton nutrient quota
"""
function net_linear_loss_quota(P, linear_mortality, DOM_POM_fractionation, quota)
    # sum over all plankton in `P`
    return sum([
        linear_loss(P[name], linear_mortality[replace(name, r"\d+" => "")]) *
        quota *
        DOM_POM_fractionation for name in names(P, 1)
    ])
end

"""
    net_quadratic_loss_quota(P, quadratic_mortality, DOM_POM_fractionation, quota, plankton_type_prefix=["Z"])

Net loss of all plankton due to quadratic mortality with a quota term.
The quota term is multiplied with estimated loss to convert from e.g. carbon to nitrogen.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `quadratic_mortality`: NamedArray of all plankton quadratic mortality rates
- `DOM_POM_fractionation`: Float which represents dissolved-to-particulate fractionation
- `quota`: Float which represents plankton nutrient quota
- `plankton_type_prefix`: Array of prefixes used in plankton names to indicate their type,
    use here to sum over only the relevant plankton (e.g., "Z" for zooplankton)
"""
function net_quadratic_loss_quota(
    P, quadratic_mortality, DOM_POM_fractionation, quota, plankton_type_prefix=["Z"]
)
    # sum over plankton that have a `quadratic_mortality`
    return sum([
        quadratic_loss(P[name], quadratic_mortality[replace(name, r"\d+" => "")]) *
        quota *
        DOM_POM_fractionation for
        name in names(P, 1) if any(prefix -> occursin(prefix, name), plankton_type_prefix)
    ])
end

end # module
