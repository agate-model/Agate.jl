using Agate.Library.Growth
using Agate.Library.Mortality
using Agate.Library.Nutrients
using Agate.Library.Photosynthesis
using Agate.Library.Predation
using Agate.Library.Remineralization

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
    prey_name, P, maximum_predation_rate, holling_half_saturation, palatability
)
    # get predator names from maximum_predation_rate array (prey has none)
    loss = sum(
        preferential_predation_loss(
            P[prey_name],
            P[predator_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[predator_name],
            palatability[predator_name, prey_name],
        ) for predator_name in names(maximum_predation_rate)[1]
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
        preferential_predation_gain(
            P[prey_name],
            P[predator_name],
            assimilation_efficiency[predator_name, prey_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[predator_name],
            palatability[predator_name, prey_name],
        ) for prey_name in names(P)[1]
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
        preferential_predation_assimilation_loss(
            P[prey_name],
            P[predator_name],
            assimilation_efficiency[predator_name, prey_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[predator_name],
            palatability[predator_name, prey_name],
        ) for prey_name in names(P)[1]
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
    return sum([linear_loss(P[name], linear_mortality[name]) for name in names(P)[1]]) * fraction
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
        [
            quadratic_loss(P[name], quadratic_mortality[name]) for
            name in names(quadratic_mortality)[1]
        ] * fraction,
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
    N, P, PAR, maximum_growth_rate, nitrogen_half_saturation, alpha
)
    return sum([
        # only phytoplankton have maximum_growth_rate, nitrogen_half_saturation and alpha
        # --> get names from either of those arrays
        idealized_photosynthetic_growth(
            N,
            P[name],
            PAR,
            maximum_growth_rate[name],
            nitrogen_half_saturation[name],
            alpha[name],
        ) for name in names(maximum_growth_rate)[1]
    ],)
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
        ) for predator_name in names(maximum_predation_rate)[1]
    ])
end
