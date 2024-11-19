using Oceananigans.Units
using Agate

include(joinpath("..", "N2P2ZD", "tracers.jl"))

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
        "alpha" => 0.1953 / day,
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
        "alpha" => 1e-99 / day,
        "assimilation_efficiency" => 0.32,
    ),
)

emergent_parameters = compute_darwin_parameters(defined_parameters)

# simple test --> compare to parameters in N2P2ZD example

plankton_order = ["P1", "P2", "Z1", "Z2"]

for (key, params) in emergent_parameters
    # start with arrays of values
    if !(key ∈ ["assimilation_efficiency", "palatability", "volume"])
        comparison = all(parameters[Symbol(key)] .== [params[p] for p in plankton_order])
        println(key, " values are the same: ", comparison)
        # matrices of values -> compare row at a time
    elseif !(key == "volume")
        for (i, p) in enumerate(plankton_order)
            emergent_row = params[p, :]
            true_row = parameters[Symbol(key)][i, :]
            comparison = all(true_row .== [emergent_row[p] for p in plankton_order])
            println(key, " ", p, " values are the same: ", comparison)
        end
    end
end

# for simplicity define the biogeochemistry dict seperately
biogeochemistry_parameters = Dict(
    "detritus_remineralization" => 0.1213 / day,
    "feeding_export_poc_doc_fraction" => 0.5,
    "mortality_export_fraction" => 0.5,
)
created_parameters = merge(biogeochemistry_parameters, emergent_parameters)
#note that this dictionary would need to be converted to a named tuple to work with create_bgc_struc()...
