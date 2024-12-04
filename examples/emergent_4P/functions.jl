# phytoplankton growth
function menden_limitation(N, nitrogen_half_saturation)
    return N / (nitrogen_half_saturation + N)
end
function smith_light_limitation(PAR, alpha, maximum_growth_rate)
    if alpha == 0 && maximum_growth_rate == 0
        return 0.0
    end
    return alpha * PAR / sqrt(maximum_growth_rate^2 + alpha^2 * PAR^2)
end
function photosynthetic_growth(
    N,
    P,
    PAR,
    maximum_growth_rate,
    nitrogen_half_saturation,
    alpha,
)
    return maximum_growth_rate *
           menden_limitation(N, nitrogen_half_saturation) *
           smith_light_limitation(PAR, alpha, maximum_growth_rate) *
           P
end
# zooplankton growth
holling_type_2(R, k) = R / (k + R)

"""
Estimates the loss rate of P (prey) to Z (predator).

# Arguments
- `P`: prey plankton
- `Z`: predator plankton
- `maximum_predation_rate`: maximum predation rate of the predator
- `holling_half_saturation`: holling half saturation constant of the predator
- `palatability`: palatability of prey to the predator
"""
function predation_loss(
    P,
    Z,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    loss =
        maximum_predation_rate *
        palatability *
        holling_type_2(P, holling_half_saturation) *
        Z
    return loss
end

"""
Estimates the gain rate of Z (predator) feeding on P (prey).

# Arguments
- `P`: prey plankton
- `Z`: predator plankton
- `maximum_predation_rate`: maximum predation rate of the predator
- `holling_half_saturation`: holling half saturation constant of the predator
- `palatability`: palatability of prey to the predator
- `assimilation_efficiency`: assimilation efficiency of prey to the predator
"""
function predation_gain(
    P,
    Z,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    gain =
        predation_loss(
            P, Z, maximum_predation_rate, holling_half_saturation, palatability
        ) * assimilation_efficiency
    return gain
end

"""
Estimates the rate at which plankton predation gain is lost due to inefficient assimilation
efficiency (e.g. 'sloppy feeding').

Usually, this is used to estimate fluxes from predation to dissolved and particulate
organic matter (DOM and POM).

# Arguments
- `P`: prey plankton
- `Z`: predator plankton
- `maximum_predation_rate`: maximum predation rate of the predator
- `holling_half_saturation`: holling half saturation constant of the predator
- `palatability`: palatability of prey to the predator
- `assimilation_efficiency`: assimilation efficiency of prey to the predator
"""
function predation_assimilation_loss(
    P,
    Z,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    assimilation_loss =
        predation_loss(
            P, Z, maximum_predation_rate, holling_half_saturation, palatability
        ) * (1 - assimilation_efficiency)
    return assimilation_loss
end

"""
Estimates the total loss rate of the prey `P[prey_name]` to predation.

For plankton `P[prey_name]`, the function loops over each predator to
estimate the total loss of plankton `prey_name` due to predation.

# Arguments
- `prey_name`: name of the prey plankton to access value as `P[prey_name]`
- `P`: NamedArray which includes all plankton
- `maximum_predation_rate`: NamedArray of all plankton predation rates
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function summed_predation_loss(
    prey_name,
    P,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    # get predator names from maximum_predation_rate array (prey has none)
    loss = sum(
        predation_loss(
            P[prey_name],
            P[predator_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[predator_name],
            palatability[predator_name, prey_name],
        ) for predator_name in names(maximum_predation_rate)
    )

    return loss
end

"""
Estimates the total predation gain of the predator (P[predator_name]) feeding on all plankton.

For plankton P[predator_name], the function loops over each prey (P[prey_name]) to
estimate the total gain due to predation.

# Arguments
- `predator_name`: name of the predator, e.g. P[predator_name]
- `P`: NamedArray which includes all plankton
- `maximum_predation_rate`: NamedArray of all plankton predation rates
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `assimilation efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function summed_predation_gain(
    predator_name,
    P,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    # sum over all plankton in P (return 0 if not suitable prey for this predator)
    gain = sum(
        predation_gain(
            P[prey_name],
            P[predator_name],
            assimilation_efficiency[predator_name, prey_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[predator_name],
            palatability[predator_name, prey_name],
        ) for prey_name in names(P)
    )

    return gain
end

"""
Estimates the total assimilation loss of the predator (P[predator_name]) feeding on all plankton.

For plankton P[predator_name], the function loops over each prey (P[prey_name]) to
estimate the total assimilation loss during predation.

# Arguments
- `predator_name`: name of the predator, e.g. P[predator_name]
- `P`: NamedArray which includes all plankton
- `maximum_predation_rate`: NamedArray of all plankton predation rates
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `assimilation efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function summed_predation_assimilation_loss(
    predator_name,
    P,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    # sum over all plankton in P (return 0 if not suitable prey for this predator)
    assimilation_loss = sum(
        predation_assimilation_loss(
            P[prey_name],
            P[predator_name],
            assimilation_efficiency[predator_name, prey_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[predator_name],
            palatability[predator_name, prey_name],
        ) for prey_name in names(P)
    )

    return assimilation_loss
end

#mortality
linear_loss(P, linear_mortality) = linear_mortality * P
quadratic_loss(P, quadratic_mortality) = quadratic_mortality * P^2
#detritus
function remineralization(D, detritus_remineralization)
    return D * detritus_remineralization
end
#sums
"""
Net loss of all plankton due to linear mortality.

# Arguments
- `P`: NamedArray which includes all plankton
- `linear_mortality`: NamedArray of all plankton linear mortality rates
"""
function net_linear_loss(P, linear_mortality, fraction)
    # all plankton have a linear loss value --> get all P names
    return sum([linear_loss(P[name], linear_mortality[name]) for name in names(P)]) * fraction
end
"""
Net loss of all plankton due to quadratic mortality.

# Arguments
- `P`: NamedArray which includes all plankton
- `quadratic_mortality`: NamedArray of all plankton quadratic mortality rates
"""
function net_quadratic_loss(P, quadratic_mortality, fraction)
    # only zooplankton have quadratic mortality --> get names from associated array
    return sum(
        [quadratic_loss(P[name], quadratic_mortality[name]) for name in names(quadratic_mortality)] * fraction
    )
end
"""
Net photosynthetic growth of all plankton.

# Arguments
- `N`: Nitrogen
- `P`: NamedArray which includes all plankton
- `PAR`: PAR
- `maximum_growth_rate`: NamedArray of all plankton maximum growth rates
- `nitrogen_half_saturation`: NamedArray of all plankton nitrogen half saturation constants
"""
function net_photosynthetic_growth(
    N,
    P,
    PAR,
    maximum_growth_rate,
    nitrogen_half_saturation,
    alpha,
)
    return sum([
        # only phytoplankton have maximum_growth_rate, nitrogen_half_saturation and alpha
        # --> get names from either of those arrays
        photosynthetic_growth(
            N, P[name], PAR, maximum_growth_rate[name], nitrogen_half_saturation[name], alpha[name]
        ) for name in names(maximum_growth_rate)
    ])
end

"""
Net predator assimilation loss of all plankton.

# Arguments
- `P`: NamedArray which includes all plankton
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `maximum_predation_rate`: NamedArray of all plankton maximum predation rates
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `assimilation efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function net_predation_assimilation_loss(
    P,
    holling_half_saturation,
    maximum_predation_rate,
    assimilation_efficiency,
    palatability,
)
    # get predator names from maximum_predation_rate array (prey has none)
    return sum([
        summed_predation_assimilation_loss(
            predator_name,
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
        ) for predator_name in names(maximum_predation_rate)
    ])
end

"""
Wrapper function to estimate the rate at which phytoplankton biomass changes over time.

# Arguments
- `plankton_name`: The name of the plankton for which the rate of change is estimated
- `N`: Nitrogen
- `P`: NamedArray which includes all plankton
- `linear_mortality`: NamedArray of all plankton linear mortality rates
- `quadratic_mortality`: ....
- `maximum_growth_rate`: NamedArray of all plankton maximum growth rates
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `maximum_predation_rate`: NamedArray of all plankton maximum predation rates
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function phytoplankton_dt(
    plankton_name,
    N,
    P,
    PAR,
    linear_mortality,
    quadratic_mortality,
    maximum_growth_rate,
    holling_half_saturation,
    nitrogen_half_saturation,
    alpha,
    maximum_predation_rate,
    palatability,
)
    growth =
        photosynthetic_growth(
            N,
            P[plankton_name],
            PAR,
            maximum_growth_rate[plankton_name],
            nitrogen_half_saturation[plankton_name],
            alpha[plankton_name],
        ) - summed_predation_loss(
            plankton_name, P, maximum_predation_rate, holling_half_saturation, palatability
        ) - linear_loss(P[plankton_name], linear_mortality[plankton_name]) -
        quadratic_loss(P[plankton_name], quadratic_mortality[plankton_name])
    return growth
end

"""
Wrapper function to estimate the rate at which zooplankton biomass changes over time.

# Arguments
- `plankton_name`: The name of the plankton for which the rate of change is estimated
- `P`: NamedArray which includes all plankton
- `linear_mortality`: NamedArray of all plankton linear mortality rates
- `quadratic_mortality`: ....
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `maximum_predation_rate`: NamedArray of all plankton maximum predation rates
- `assimilation efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function zooplankton_dt(
    plankton_name,
    P,
    linear_mortality,
    quadratic_mortality,
    holling_half_saturation,
    maximum_predation_rate,
    assimilation_efficiency,
    palatability,
)
    growth = summed_predation_gain(
            plankton_name,
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
        ) - linear_loss(P[plankton_name], linear_mortality[plankton_name]) -
        quadratic_loss(P[plankton_name], quadratic_mortality[plankton_name])
    return growth
end
