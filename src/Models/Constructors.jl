module Constructors

include("Biogeochemistry.jl")
include("Parameters.jl")
include("Tracers.jl")

using .Biogeochemistry
using .Parameters
using .Tracers

using Oceananigans.Units

export construct_NPZD_instance

function construct_NPZD_instance(
    n_phyto=2,
    n_zoo=2,
    nutrients_dt=typical_nutrients,
    detritus_dt=typical_detritus,
    phyto_growth=simplified_phytoplankton_growth,
    zoo_growth=simplified_zooplankton_growth,
    phyto_args=Dict(
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
    zoo_args=Dict(
        "volumes" => Dict(
            "min_volume" => 10, "max_volume" => 100, "splitting" => "linear_splitting"
        ),
        "allometry" => Dict(
            "maximum_growth_rate" => Dict("a" => 0, "b" => 0),
            "nitrogen_half_saturation" => Dict("a" => 0, "b" => 0),
            "maximum_predation_rate" => Dict("a" => 1, "b" => 1),
        ),
        "palatability" => Dict("optimum_predator_prey_ratio" => 10, "protection" => 1),
        "assimilation_efficiency" =>
            Dict("can_be_eaten" => 0, "can_eat" => 1, "assimilation_efficiency" => 0.32),
        # anything else goes here
        "linear_mortality" => 8e-7 / second,
        "holling_half_saturation" => 5.0,
        "quadratic_mortality" => 1e-6 / second,
        "alpha" => 0,
    ),
    bgc_args=Dict(
        "detritus_remineralization" => 0.1213 / day,
        "feeding_export_poc_doc_fraction" => 0.5,
        "mortality_export_fraction" => 0.5,
    ),
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
)
    phyto_args["n"] = n_phyto
    zoo_args["n"] = n_zoo

    # TODO: if palatability_matrix or assimilation_efficiency_matrix are defined
    # (not nothing) remove them from the dictionaries

    defined_parameters = Dict("P" => phyto_args, "Z" => zoo_args)
    emergent_parameters = compute_darwin_parameters(defined_parameters)
    parameters = NamedTuple(
        Symbol(k) => v for (k, v) in merge(bgc_args, emergent_parameters)
    )

    plankton_array = vcat(
        [Symbol("P$i") for i in 1:n_phyto], [Symbol("Z$i") for i in 1:n_zoo]
    )
    tracers = Dict("N" => nutrients_dt(plankton_array), "D" => detritus_dt(plankton_array))
    for i in 1:n_phyto
        name = "P$i"
        tracers[name] = phyto_growth(plankton_array, name)
    end
    for i in 1:n_zoo
        name = "Z$i"
        tracers[name] = zoo_growth(plankton_array, name)
    end

    return define_tracer_functions(parameters, tracers)
end

end # module
