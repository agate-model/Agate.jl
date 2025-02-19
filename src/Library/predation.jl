"""
Functions related to zooplankton predation.
"""

module Predation

export holling_type_2,
    predation_loss_idealized,
    predation_gain_idealized,
    predation_assimilation_loss_idealized,
    predation_loss_preferential,
    predation_gain_preferential,
    predation_assimilation_loss_preferential,
    summed_predation_gain_preferential,
    summed_predation_loss_preferential,
    summed_predation_assimilation_loss_preferential,
    net_predation_assimilation_loss_preferential,
    net_predation_assimilation_loss_preferential_fractionated,
    net_predation_assimilation_loss_preferential_fractionated_quota,
    assimilation_efficiency_emergent_binary

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
function predation_loss_idealized(P, Z, gₘₐₓ, kₚ)
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
function predation_gain_idealized(P, Z, β, gₘₐₓ, kₚ)
    return β * predation_loss_idealized(P, Z, gₘₐₓ, kₚ)
end

"""
Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

Note that this differs from the predation_gain_idealized as this function represents the transfer of biomass from the prey to the environment
rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing/holling half saturation
"""
function predation_assimilation_loss_idealized(P, Z, β, gₘₐₓ, kₚ)
    return (1 - β) * predation_loss_idealized(P, Z, gₘₐₓ, kₚ)
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
function predation_loss_preferential(P, Z, gₘₐₓ, kₚ, palatability)
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
function predation_gain_preferential(P, Z, β, gₘₐₓ, kₚ, palatability)
    return β * predation_loss_preferential(P, Z, gₘₐₓ, kₚ, palatability)
end

"""
Estimates the rate at which plankton predation gain is lost to the environment due to inefficient assimilation efficiency
(e.g. 'sloppy feeding').

Note that this differs from the predation_gain_preferential as this function represents the transfer of biomass from the prey to the environment
rather than the transfer of biomass from the prey to the predator.

# Arguments
- `P`: phytoplankton concentration
- `Z`: zooplankton concentration
- `β`: assimilation efficiency of prey to the predator
- `gₘₐₓ`: maximum grazing rate of the predator
- `kₚ`: grazing/holling half saturation
- `palatability`: the likelihood at which the predator feeds on the prey
"""
function predation_assimilation_loss_preferential(P, Z, β, gₘₐₓ, kₚ, palatability)
    return (1 - β) * predation_loss_preferential(P, Z, gₘₐₓ, kₚ, palatability)
end

"""
Estimates the total loss rate of the prey `P[prey_name]` to predation.

For plankton `P[prey_name]`, the function loops over each predator to
estimate the total loss of plankton `prey_name` due to predation.

# Arguments
- `prey_name`: name of the prey plankton to access value as `P[prey_name]`
- `P`: NamedArray which includes all plankton
- `maximum_predation_rate`: NamedArray of all plankton predation rates
- `holling_half_saturation`: predation half saturation constant
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `plankton_type_prefix`: Array of prefixes used in plankton names to indicate their type,
    use here to sum over only predator plankton (e.g., "Z" for zooplankton)
"""
function summed_predation_loss_preferential(
    prey_name,
    P,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
    plankton_type_prefix=["Z"],
)
    loss = sum(
        predation_loss_preferential(
            P[prey_name],
            P[predator_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[replace(predator_name, r"\d+" => "")],
            palatability[predator_name, prey_name],
        ) for predator_name in names(P, 1) if
        any(prefix -> occursin(prefix, predator_name), plankton_type_prefix)
    )

    return loss
end

"""
Estimates the total predation gain of the predator (`P[predator_name]`) feeding on all plankton.

For plankton `P[predator_name]`, the function loops over each prey (`P[prey_name]`) to
estimate the total gain due to predation.

# Arguments
- `predator_name`: name of the predator, e.g. `P[predator_name]`
- `P`: NamedArray which includes all plankton concentration values
- `maximum_predation_rate`: NamedArray of all plankton predation rates
- `holling_half_saturation`: predation half saturation constant
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `assimilation_efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function summed_predation_gain_preferential(
    predator_name,
    P,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    # sum over all plankton in P (return 0 if not suitable prey for this predator)
    gain = sum(
        predation_gain_preferential(
            P[prey_name],
            P[predator_name],
            assimilation_efficiency[predator_name, prey_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[replace(predator_name, r"\d+" => "")],
            palatability[predator_name, prey_name],
        ) for prey_name in names(P, 1)
    )

    return gain
end

"""
Estimates the total assimilation loss of the predator (`P[predator_name]`) feeding on all plankton.

For plankton P`[predator_name]`, the function loops over each prey (`P[prey_name]`) to
estimate the total assimilation loss during predation.

# Arguments
- `predator_name`: name of the predator, e.g. `P[predator_name]`
- `P`: NamedArray which includes all plankton concentration values
- `maximum_predation_rate`: NamedArray of all plankton predation rates
- `holling_half_saturation`: plankton predation half saturation constant
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `assimilation_efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function summed_predation_assimilation_loss_preferential(
    predator_name,
    P,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
)
    # sum over all plankton in P (return 0 if not suitable prey for this predator)
    assimilation_loss = sum(
        predation_assimilation_loss_preferential(
            P[prey_name],
            P[predator_name],
            assimilation_efficiency[predator_name, prey_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[replace(predator_name, r"\d+" => "")],
            palatability[predator_name, prey_name],
        ) for prey_name in names(P, 1)
    )

    return assimilation_loss
end

"""
Estimates the total assimilation loss of the predator (`P[predator_name]`) feeding on all plankton.

For plankton P`[predator_name]`, the function loops over each prey (`P[prey_name]`) to
estimate the total assimilation loss during predation.

# Arguments
- `predator_name`: name of the predator, e.g. `P[predator_name]`
- `P`: NamedArray which includes all plankton concentration values
- `maximum_predation_rate`: NamedArray of all plankton predation rates
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `assimilation_efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
"""
function summed_predation_assimilation_loss_preferential_quota(
    predator_name,
    P,
    assimilation_efficiency,
    maximum_predation_rate,
    holling_half_saturation,
    palatability,
    quota,
)
    # sum over all plankton in P (return 0 if not suitable prey for this predator)
    assimilation_loss = sum(
        predation_assimilation_loss_preferential(
            P[prey_name],
            P[predator_name],
            assimilation_efficiency[predator_name, prey_name],
            maximum_predation_rate[predator_name],
            holling_half_saturation[replace(predator_name, r"\d+" => "")],
            palatability[predator_name, prey_name],
        ) * quota[replace(prey_name, r"\d+" => "")] for prey_name in names(P, 1)
    )

    return assimilation_loss
end

"""
Net predator assimilation loss of all plankton.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `maximum_predation_rate`: NamedArray of all plankton maximum predation rates
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `assimilation_efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `plankton_type_prefix`: Array of prefixes used in plankton names to indicate their type,
    use here to sum over only predator plankton (e.g., "Z" for zooplankton)
"""
function net_predation_assimilation_loss_preferential(
    P,
    holling_half_saturation,
    maximum_predation_rate,
    assimilation_efficiency,
    palatability,
    plankton_type_prefix=["Z"],
)
    # get predator names from `maximum_predation_rate` array (prey has none)
    return sum([
        summed_predation_assimilation_loss_preferential(
            predator_name,
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
        ) for predator_name in names(P, 1) if
        any(prefix -> occursin(prefix, predator_name), plankton_type_prefix)
    ])
end

"""
Net predator assimilation loss of all plankton which is fractionated between DOC and POC.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `maximum_predation_rate`: NamedArray of all plankton maximum predation rates
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `assimilation_efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `DOM_POM_fractionation`: float representing the fraction of loss going to DOM and POM.
"""
function net_predation_assimilation_loss_preferential_fractionated(
    P,
    holling_half_saturation,
    maximum_predation_rate,
    assimilation_efficiency,
    palatability,
    DOM_POM_fractionation,
    plankton_type_prefix=["Z"],
)
    # get predator names from `maximum_predation_rate` array (prey has none)
    return sum([
        summed_predation_assimilation_loss_preferential(
            predator_name,
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
        ) for predator_name in names(maximum_predation_rate, 1) if
        any(prefix -> occursin(prefix, predator_name), plankton_type_prefix)
    ]) * DOM_POM_fractionation
end

"""
Net predator assimilation loss of all plankton which is fractionated between DOC and POC.
With a quota term.

# Arguments
- `P`: NamedArray which includes all plankton concentration values
- `holling_half_saturation`: NamedArray of all plankton predation half saturation constants
- `maximum_predation_rate`: NamedArray of all plankton maximum predation rates
- `palatability`: NamedArray of all plankton palatabilities where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `assimilation_efficiency`: NamedArray of all plankton assimilation efficiencies where:
    - each row is a predator
    - each column is a prey
    - values are accessed as `palat[predator, prey]`
    - for a non-predator [i,:]=0
- `DOM_POM_fractionation`: float representing the fraction of loss going to DOM and POM.
- `quota`: NamedArray of all plankton quotas


"""
function net_predation_assimilation_loss_preferential_fractionated_quota(
    P,
    holling_half_saturation,
    maximum_predation_rate,
    assimilation_efficiency,
    palatability,
    DOM_POM_fractionation,
    quota,
    plankton_type_prefix=["Z"],
)
    # get predator names from `maximum_predation_rate` array (prey has none)
    return sum([
        summed_predation_assimilation_loss_preferential_quota(
            predator_name,
            P,
            assimilation_efficiency,
            maximum_predation_rate,
            holling_half_saturation,
            palatability,
            quota,
        ) for predator_name in names(maximum_predation_rate, 1) if
        any(prefix -> occursin(prefix, predator_name), plankton_type_prefix)
    ]) * DOM_POM_fractionation
end

"""
    assimilation_efficiency_emergent_binary(prey_data, predator_data)

Determines the assimilation efficiency of a predator consuming prey, based on binary conditions of edibility.

The function evaluates whether the predator can eat the prey and whether the prey can be consumed, and assigns the assimilation efficiency accordingly.

# Arguments
- `prey_data`: A dictionary containing prey-specific data:
  - `can_be_eaten`: A binary value (1 or 0) indicating if the prey can be consumed by the predator.
- `predator_data`: A dictionary containing predator-specific data:
  - `can_eat`: A binary value (1 or 0) indicating if the predator can consume prey.
  - `assimilation_efficiency`: The efficiency with which the predator assimilates nutrients from the prey if the conditions are met.

# Returns
- `assimilation_efficiency`:
  - If `can_eat` is 1 and `can_be_eaten` is 1, returns the predator's `assimilation_efficiency`.
  - Otherwise, returns 0.
"""
function assimilation_efficiency_emergent_binary(prey_data, predator_data)
    if predator_data["can_eat"] == 1 && prey_data["can_be_eaten"] == 1
        assimilation_efficiency = predator_data["assimilation_efficiency"]
    elseif predator_data["can_eat"] == 1 && prey_data["can_be_eaten"] == 0
        assimilation_efficiency = 0
    elseif predator_data["can_eat"] == 0
        assimilation_efficiency = 0
    end
    return assimilation_efficiency
end

end # module
