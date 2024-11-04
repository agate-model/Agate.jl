include("emergent_functions.jl")
include("helper_functions.jl")
# Sample plankton data with required keys
plankton = Dict(
    "species1" => Dict(
        "volume" => 1.0,
        "volume_a" => 1.0,
        "volume_b" => 1.0,
        "pred" => 2.0,
        "optimum_predator_prey_ratio" => 10,
        "protection" => 0,
        "growth_a" => 1,
        "growth_b" => 1,
    ),
    "species2" => Dict(
        "volume" => 1.5,
        "volume_a" => 1.0,
        "volume_b" => 1.0,
        "pred" => 1.8,
        "optimum_predator_prey_ratio" => 10,
        "protection" => 0,
        "growth_a" => 1,
        "growth_b" => 1,
    ),
)

# Dictionary of emergent functions with symbolic expressions
emergent_functions = Dict(
    "growth_rate" => (dummy_emergent_growth, ["growth_a", "growth_b", "volume"]),
    "palatability" =>
        (dummy_emergent_palat, ["volume", "optimum_predator_prey_ratio", "protection"]),
    "predation_rate" => (dummy_emergent_predation_rate, ["volume_a", "volume_b", "volume"]),
    "nitrogen_half_saturation" =>
        (dummy_emergent_nitrogen_half_saturation, ["volume_a", "volume_b", "volume"]),
)

# Dictionary to store results
results_dict = Dict()

# Main loop: analyze each function in `emergent_functions` dynamically and store the results
for (func_name, (func, param_names)) in emergent_functions
    println("Analyzing function: $func_name")

    # Run the analysis and get the result
    result = emergent_analysis(plankton, func, param_names)

    # Store the result in the dictionary
    results_dict[func_name] = result
end

# Display the results dictionary
println("All results:", results_dict)
