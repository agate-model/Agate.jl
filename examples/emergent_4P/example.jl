using Oceananigans.Units
using Agate

defined_parameters = Dict(
    "P" => Dict(
        "n" => 2,
        "min_volume" => 1,
        "max_volume" => 10,
        "splitting" => "log_splitting",
        "growth_a" => 1,
        "growth_b" => 1,
        "protection" => 0,
        "optimum_predator_prey_ratio" => 0,
        "nitrogen_half_saturation_a" => 1,
        "nitrogen_half_saturation_b" => 1,
        "predation_rate_a" => 0,
        "predation_rate_b" => 0,
        "linear_mortality" => 8e-7 / second,
        "holling_half_saturation" => 0,
        "quadratic_mortality" => 0,
        "alpha" => 0.1953,
        "assimilation_efficiency" => 0,
    ),
    "Z" => Dict(
        "n" => 2,
        "min_volume" => 10,
        "max_volume" => 100,
        "splitting" => "linear_splitting",
        "growth_a" => 0,
        "growth_b" => 0,
        "protection" => 1,
        "optimum_predator_prey_ratio" => 10,
        "nitrogen_half_saturation_a" => 0,
        "nitrogen_half_saturation_b" => 0,
        "predation_rate_a" => 1,
        "predation_rate_b" => 1,
        "linear_mortality" => 8e-7 / second,
        "holling_half_saturation" => 5.0,
        "quadratic_mortality" => 1e-6 / second,
        "alpha" => 1e-99,
        "assimilation_efficiency" => 0.32,
    ),
)

emergent_parameters = compute_darwin_parameters(defined_parameters)

# for simplicity define the biogeochemistry dict seperately
biogeochemistry_parameters = Dict(
    "detritus_remineralization" => 0.1213 / day,
    "feeding_export_poc_doc_fraction" => 0.5,
    "mortality_export_fraction" => 0.5,
)

#merge into one dictionary
parameters = merge(biogeochemistry_parameters, emergent_parameters)

println(parameters)

# compare to parameters in N2P2ZD example

#note that this dictionary would need to be converted to a named tuple to work with create_bgc_struc()...
