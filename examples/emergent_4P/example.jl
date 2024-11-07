include("emergent_functions.jl")
include("helper_functions.jl")

intermediate_parameters = Dict(
    "P1" => Dict(
        "volume" => 1,
        "growth_a" => 1,
        "growth_b" => 1,
        "protection" => 0,
        "optimum_predator_prey_ratio" => 0,
        "nitrogen_half_saturation_a" => 1,
        "nitrogen_half_saturation_b" => 1,
        "predation_rate_a" => 0,
        "predation_rate_b" => 0,
    ),
    "P2" => Dict(
        "volume" => 10,
        "growth_a" => 1,
        "growth_b" => 1,
        "protection" => 0,
        "optimum_predator_prey_ratio" => 0,
        "nitrogen_half_saturation_a" => 1,
        "nitrogen_half_saturation_b" => 1,
        "predation_rate_a" => 0,
        "predation_rate_b" => 0,
    ),
    "Z1" => Dict(
        "volume" => 10,
        "growth_a" => 0,
        "growth_b" => 0,
        "protection" => 1,
        "optimum_predator_prey_ratio" => 10,
        "nitrogen_half_saturation_a" => 0,
        "nitrogen_half_saturation_b" => 0,
        "predation_rate_a" => 1,
        "predation_rate_b" => 1,
    ),
    "Z2" => Dict(
        "volume" => 100,
        "growth_a" => 0,
        "growth_b" => 0,
        "protection" => 1,
        "optimum_predator_prey_ratio" => 10,
        "nitrogen_half_saturation_a" => 0,
        "nitrogen_half_saturation_b" => 0,
        "predation_rate_a" => 1,
        "predation_rate_b" => 1,
    ),
)

# Dictionary of emergent functions with symbolic expressions
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

# Dictionary to store results
emergent_parameters = Dict()

# Main loop: analyze each function in `emergent_functions` dynamically and store the results
for (func_name, (func, param_names)) in emergent_functions
    println("Analyzing function: $func_name")

    # Run the analysis and get the result
    result = emergent_analysis(intermediate_parameters, func, param_names)

    # Store the result in the dictionary
    emergent_parameters[func_name] = result
end

# Display the results dictionary
println("All results:", emergent_parameters)

#check named arrays:

#growth rate of P2:
println(emergent_parameters["growth_rate"]["P2"])

#palability of P2 to P1:
println(emergent_parameters["palatability"]["P1", "P2"])

#palability of Z1 to P1:
println(emergent_parameters["palatability"]["Z1", "P2"])
