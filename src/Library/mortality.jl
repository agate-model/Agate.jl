module Mortality

export linear_loss, quadratic_loss, net_linear_loss, net_quadratic_loss

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
Net loss of all plankton due to linear mortality.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `linear_mortality`: NamedArray of plankton linear mortality rates
"""
function net_linear_loss(P, linear_mortality, fraction)
    # sum over all plankton in `P` - strip digits from plankton name to get its type (e.g., "P")
    return sum([
        linear_loss(P[name], linear_mortality[replace(name, r"\d+" => "")]) for
        name in names(P, 1)
    ]) * fraction
end

"""
Net loss of all plankton due to quadratic mortality.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `quadratic_mortality`: plankton quadratic mortality rate
- `plankton_type_prefix`: Array of prefixes used in plankton names to indicate their type,
    use here to sum over only the relevant plankton (e.g., "Z" for zooplankton)
"""
function net_quadratic_loss(P, quadratic_mortality, fraction, plankton_type_prefix=["Z"])
    return sum(
        [
            quadratic_loss(P[name], quadratic_mortality) for name in names(P, 1) if
            any(prefix -> occursin(prefix, name), plankton_type_prefix)
        ] * fraction,
    )
end

end # module
