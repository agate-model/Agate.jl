using Oceananigans.Units

include("emergent_functions.jl")
include("helper_functions.jl")

defined_parameters = Dict(
    "P" => Dict(
        "n" => 2,
        "min_volume" => 1,
        "max_volume" => 10,
        "splitting" => log_splitting,
        "growth_a" => 1,
        "growth_b" => 1,
        "protection" => 0,
        "optimum_predator_prey_ratio" => 0,
        "nitrogen_half_saturation_a" => 1,
        "nitrogen_half_saturation_b" => 1,
        "predation_rate_a" => 0,
        "predation_rate_b" => 0,
        "linear_mortality" => 8e-7 / second,
    ),
    "Z" => Dict(
        "n" => 2,
        "min_volume" => 10,
        "max_volume" => 100,
        "splitting" => linear_splitting,
        "growth_a" => 0,
        "growth_b" => 0,
        "protection" => 1,
        "optimum_predator_prey_ratio" => 10,
        "nitrogen_half_saturation_a" => 0,
        "nitrogen_half_saturation_b" => 0,
        "predation_rate_a" => 1,
        "predation_rate_b" => 1,
        "linear_mortality" => 8e-7 / second,
    ),
)

# Generate intermediate parameters 
intermediate_parameters = split_size_parameters(defined_parameters)

# Dictionary of emergent functions
emergent_functions = Dict(
    "growth_rate" => (dummy_emergent_growth, ["growth_a", "growth_b", "volume"]),
    "palatability" =>
        (dummy_emergent_palat, ["volume", "optimum_predator_prey_ratio", "protection"]),
    "predation_rate" =>
        (dummy_emergent_predation_rate, ["predation_rate_a", "predation_rate_b", "volume"]),
    "nitrogen_half_saturation" => (
        dummy_emergent_nitrogen_half_saturation,
        ["nitrogen_half_saturation_a", "nitrogen_half_saturation_b", "volume"],
    ),
)

emergent_parameters = analyze_and_merge(emergent_functions, intermediate_parameters)

# Display the results dictionary
println("All results:", emergent_parameters)

#check named arrays:

#growth rate of P2:
println(emergent_parameters["growth_rate"]["P2"])

#palability of P2 to P1:
println(emergent_parameters["palatability"]["P1", "P2"])

#palability of Z1 to P1:
println(emergent_parameters["palatability"]["Z1", "P2"])

#for simplicity define the biogeochemistry dict seperately

biogeochemistry_parameters = Dict(
    "detritus_remineralization" => 0.1213 / day,
    "feeding_export_poc_doc_fraction" => 0.5,
    "mortality_export_fraction" => 0.5,
)

parameters = merge(biogeochemistry_parameters, emergent_parameters)

print(parameters)
