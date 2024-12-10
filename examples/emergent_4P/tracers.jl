using Agate
using Agate.Models.Tracers

using NamedArrays
using Oceananigans.Units

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
        "alpha" => 0,
    ),
)

plankton_array = [:P1, :P2, :Z1, :Z2]
tracers = Dict(
    "N" => typical_nutrients(plankton_array),
    "D" => typical_detritus(plankton_array),
    "P1" => simplified_phytoplankton_growth(plankton_array, "P1"),
    "P2" => simplified_phytoplankton_growth(plankton_array, "P2"),
    "Z1" => simplified_zooplankton_growth(plankton_array, "Z1"),
    "Z2" => simplified_zooplankton_growth(plankton_array, "Z2"),
)

# merge everything and convert to NamedTuple
emergent_parameters = compute_darwin_parameters(defined_parameters)

biogeochemistry_parameters = Dict(
    "detritus_remineralization" => 0.1213 / day,
    "feeding_export_poc_doc_fraction" => 0.5,
    "mortality_export_fraction" => 0.5,
)
parameters = NamedTuple(
    Symbol(k) => v for (k, v) in merge(biogeochemistry_parameters, emergent_parameters)
)
N2P2ZD = define_tracer_functions(parameters, tracers; helper_functions="functions.jl")
