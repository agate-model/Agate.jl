using Oceananigans.Units
using Agate

# Q: SHOULD THIS BE A NAMEDTUPLE INSTEAD?
# TODO: users shouldn't have to specify 0 parameters --> do we just assume 0 and fill it OR can we update tracers to do without them entirely?
# NOTE: assimilation efficiency has to create an array of values rather than compute anything --> is there a better way to do that?
# NOTE: don't need to declare units as per second since this is the default time unit i.e., 1 == 1/sec

defined_parameters = Dict(
    "P" => Dict(
        "n" => 2,
        "volumes" =>
            Dict("min_volume" => 1, "max_volume" => 10, "splitting" => "log_splitting"),
        "allometry" => Dict(
            "maximum_growth_rate" => Dict("a" => 1, "b" => 1),
            "nitrogen_half_saturation" => Dict("a" => 1, "b" => 1),
            "maximum_predation_rate" => Dict("a" => 0, "b" => 0),
        ),
        "palatability" => Dict("optimum_predator_prey_ratio" => 0, "protection" => 0),
        "assimilation_efficiency" =>
            Dict("can_be_eaten" => 1, "can_eat" => 0, "assimilation_efficiency" => 0),
        # anything else goes here
        "linear_mortality" => 8e-7 / second,
        "holling_half_saturation" => 0,
        "quadratic_mortality" => 0,
        "alpha" => 0.1953 / day,
    ),
    "Z" => Dict(
        "n" => 2,
        "volumes" => Dict(
            "min_volume" => 10,
            "max_volume" => 100,
            "splitting" => "linear_splitting",
        ),
        "allometry" => Dict(
            "maximum_growth_rate" => Dict("a" => 0, "b" => 0),
            "nitrogen_half_saturation" => Dict("a" => 0, "b" => 0),
            "maximum_predation_rate" => Dict("a" => 1, "b" => 1),
        ),
        "palatability" => Dict("optimum_predator_prey_ratio" => 10, "protection" => 1),
        "assimilation_efficiency" => Dict(
            "can_be_eaten" => 0, "can_eat" => 1, "assimilation_efficiency" => 0.32
        ),
        # anything else goes here
        "linear_mortality" => 8e-7 / second,
        "holling_half_saturation" => 5.0,
        "quadratic_mortality" => 1e-6 / second,
        # this should be 0
        "alpha" => 1e-99 / day,
    ),
)

emergent_parameters = compute_darwin_parameters(defined_parameters)
println(emergent_parameters)

# sanity check --> compare generated values to parameters in N2P2ZD example
# NOTE: eventually the emergent functions will be updated to true ones
# that's why we don't include this in tests

include(joinpath("..", "N2P2ZD", "tracers.jl"))

plankton_order = ["P1", "P2", "Z1", "Z2"]

for (key, params) in emergent_parameters
    # start with arrays of values
    if !(key ∈ ["assimilation_efficiency_matrix", "palatability_matrix", "volumes"])
        comparison = all(parameters[Symbol(key)] .== [params[p] for p in plankton_order])
        println(key, " values are the same: ", comparison)
        # matrices of values -> compare row at a time
    elseif !(key == "volumes")
        for (i, p) in enumerate(plankton_order)
            emergent_row = params[p, :]
            true_row = parameters[Symbol(replace(key, "_matrix" => ""))][i, :]
            comparison = all(true_row .== [emergent_row[p] for p in plankton_order])
            println(key, " ", p, " values are the same: ", comparison)
        end
    end
end

# # for simplicity define the biogeochemistry dict seperately
# biogeochemistry_parameters = Dict(
#     "detritus_remineralization" => 0.1213 / day,
#     "feeding_export_poc_doc_fraction" => 0.5,
#     "mortality_export_fraction" => 0.5,
# )
# created_parameters = merge(biogeochemistry_parameters, emergent_parameters)
#note that this dictionary would need to be converted to a named tuple to work with create_bgc_struc()...
