using Agate

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
    "N" => :(
        net_linear_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            linear_mortality,
            mortality_export_fraction,
        ) +
        net_quadratic_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            quadratic_mortality,
            mortality_export_fraction,
        ) +
        remineralization(D, detritus_remineralization) - net_photosynthetic_growth(
            N,
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            PAR,
            maximum_growth_rate,
            nitrogen_half_saturation,
            alpha,
        )
    ),
    "D" => :(
        net_linear_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            linear_mortality,
            1 - mortality_export_fraction,
        ) +
        net_predation_assimilation_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            holling_half_saturation,
            maximum_predation_rate,
            assimilation_efficiency_matrix,
            palatability_matrix,
        ) +
        net_quadratic_loss(
            NamedArray([$(plankton_array...)], $(String.(plankton_array))),
            quadratic_mortality,
            1 - mortality_export_fraction,
        ) - remineralization(D, detritus_remineralization)
    ),
    "P1" => :(phytoplankton_dt(
        "P1",
        N,
        NamedArray([$(plankton_array...)], $(String.(plankton_array))),
        PAR,
        linear_mortality,
        quadratic_mortality,
        maximum_growth_rate,
        holling_half_saturation,
        nitrogen_half_saturation,
        alpha,
        maximum_predation_rate,
        palatability_matrix,
    )),
    "P2" => :(phytoplankton_dt(
        "P2",
        N,
        NamedArray([$(plankton_array...)], $(String.(plankton_array))),
        PAR,
        linear_mortality,
        quadratic_mortality,
        maximum_growth_rate,
        holling_half_saturation,
        nitrogen_half_saturation,
        alpha,
        maximum_predation_rate,
        palatability_matrix,
    )),
    "Z1" => :(zooplankton_dt(
        "Z1",
        NamedArray([$(plankton_array...)], $(String.(plankton_array))),
        linear_mortality,
        quadratic_mortality,
        holling_half_saturation,
        maximum_predation_rate,
        assimilation_efficiency_matrix,
        palatability_matrix,
    )),
    "Z2" => :(zooplankton_dt(
        "Z2",
        NamedArray([$(plankton_array...)], $(String.(plankton_array))),
        linear_mortality,
        quadratic_mortality,
        holling_half_saturation,
        maximum_predation_rate,
        assimilation_efficiency_matrix,
        palatability_matrix,
    )),
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
