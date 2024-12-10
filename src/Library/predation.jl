"""
Functions related to zooplankton predation.
"""

module Predation

export holling_type_2,
    idealized_predation_loss,
    idealized_predation_gain,
    idealized_predation_assimilation_loss,
    preferential_predation_loss,
    preferential_predation_gain,
    preferential_predation_assimilation_loss,
    summed_preferential_predation_gain,
    summed_preferential_predation_loss,
    summed_preferential_predation_assimilation_loss,
    net_preferential_predation_assimilation_loss

"""
Holling's "type II" functional response as describe in Holling 1959.
The function is similar to the Monod equation and Michaelis-Menten equation of for enzyme kinetics.
The formulation is characterized by decelerating predation as prey concentrations increase.

# Arguments
- `R`: prey density
- `k`: prey density at which predation is half it's maximum rate

"""
function holling_type_2(R::Real, k::Real)
    return R / (k + R)
end

#idealized predation (OceanBioME NPZD model)

"""
Estimates the loss rate of P (prey), to Z (predator).
In this formulation predator-prey interactions are modulated both by their density (Holling type 2)
and the prey-predator palatability.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: prey density at which predation is half it's maximum rate
"""
function idealized_predation_loss(P, Z, gₘₐₓ, kₚ)
    return gₘₐₓ * holling_type_2(P^2, kₚ^2) * Z
end

"""
Estimates the gain rate of Z (predator) feeding on P (prey).
In this formulation predation rate is multiplied by an assimilation efficiency (β) which
represents 'sloppy feeding'.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: grazing/holling half saturation
"""
function idealized_predation_gain(P, Z, β, gₘₐₓ, kₚ)
    return β * idealized_predation_loss(P, Z, gₘₐₓ, kₚ)
end

"""
Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

Note that this differs from the idealized_predation_gain as this function represents the transfer of biomass from the prey to the environment
rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing/holling half saturation
"""
function idealized_predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ)
    return (1 - β) * idealized_predation_loss(P, Z, gₘₐₓ, kₚ)
end

#preferential predation (intermediate complexity model)

"""
Estimates the loss rate of P (prey), to Z (predator).
In this formulation predator-prey interactions are modulated both by their density (Holling type 2)
and the prey-predator palatability.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: prey density at which predation is half it's maximum rate
- `palatability`: the likelihood at which the predator feeds on the prey
"""
function preferential_predation_loss(P, Z, gₘₐₓ, kₚ, palatability)
    return gₘₐₓ * palatability * holling_type_2(P, kₚ) * Z
end

"""
Estimates the gain rate of Z (predator) feeding on P (prey).
In this formulation predation rate is multiplied by an assimilation efficiency (β) which
represents 'sloppy feeding'.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency
- `gₘₐₓ`: maximum grazing rate
- `kₚ`: grazing/holling half saturation
- `palatability`: the likelihood at which the predator feeds on the prey
"""
function preferential_predation_gain(P, Z, β, gₘₐₓ, kₚ, palatability)
    return β * preferential_predation_loss(P, Z, gₘₐₓ, kₚ, palatability)
end

"""
Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

Note that this differs from the preferential_predation_gain as this function represents the transfer of biomass from the prey to the environment
rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing/holling half saturation
- `palatability`: ...
"""
function preferential_predation_assimilation_loss(P, Z, β, gₘₐₓ, kₚ, palatability)
    return (1 - β) * preferential_predation_loss(P, Z, gₘₐₓ, kₚ, palatability)
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
function summed_preferential_predation_loss(
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
function summed_preferential_predation_gain(
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
function summed_preferential_predation_assimilation_loss(
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
function net_preferential_predation_assimilation_loss(
    P,
    holling_half_saturation,
    maximum_predation_rate,
    assimilation_efficiency,
    palatability,
)
    # get predator names from maximum_predation_rate array (prey has none)
    return sum([
        summed_preferential_predation_assimilation_loss(
            predator_name,
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
        ) for predator_name in names(maximum_predation_rate)[1]
    ])
end

end # module
