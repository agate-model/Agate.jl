module Constructors

include("Biogeochemistry.jl")
include("Parameters.jl")
include("Tracers.jl")

using .Biogeochemistry
using .Parameters
using .Tracers

using Oceananigans.Units

export construct_NPZD_instance

"""
Construct an instance of an size-structured NPZD model.

This constructor builds a size-structured plankton model with two plankton functional types: phytoplankton (P) and zooplankton (Z), each of which 
can be specified to have any number of size classes (`n_phyto` and `n_zoo`). In addition to plankton, the constructor implements idealized detritus (D) and 
nitrogen (N) cycling by default, although more complex N and D cycling can also be defined using the `nutrient_dynamics` and `detritus_dynamics` arguments.

During model construction, the size of each plankton determines photosynthetic growth rates, nitrogen half saturation constants, predation rates, and optionally 
predator-prey assimilation and palatability values. Alternatively, if manually defined predator-prey assimilation and palatability values are desired, these can 
be defined using the `palatability_matrix` and `assimilation_efficiency_matrix` arguments.

Note that if non-default `*_dynamics` expressions are passed, the relevant `*_args` also need to be specified.

# Arguments
- `n_phyto`: number of phytoplankton to include in the model
- `n_zoo`: number of zooplankton to include in the model
- `nutrient_dynamics`: expression describing how nutrients change over time
- `detritus_dynamics`: expression describing how detritus evolves over time
- `phyto_dynamics`: expression describing how phytoplankton grow
- `zoo_dynamics`: expression describing how zooplankton grow
- `phyto_args`: dictionary of Phytoplankton parameters
- `zoo_args`: Dictionary of zooplankton parameters
- `palatability_args`: Dictionary of arguments from which a palatability matrix between all
   plankton can be computed
- `assimilation_effificency_args`: Dictionary of arguments from which an assimilation
   efficiency matrix between all plankton can be computed
- `bgc_args`: biogeochemistry parameters related to nutrient and detritus
- `palatability_matrix`: optional palatability matrix passed as a NamedArray, if provided then
   `paralatability_args` are ignored
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as a
   NamedArray, if provided then `assimilation_args` are ignored
"""
function construct_NPZD_instance(
    n_phyto=2,
    n_zoo=2,
    nutrient_dynamics=typical_nutrients,
    detritus_dynamics=typical_detritus,
    phyto_dynamics=simplified_phytoplankton_growth,
    zoo_dynamics=simplified_zooplankton_growth,
    phyto_args=Dict(
        "volumes" =>
            Dict("min_volume" => 1, "max_volume" => 10, "splitting" => "log_splitting"),
        "allometry" => Dict(
            "maximum_growth_rate" => Dict("a" => 1, "b" => 1),
            "nitrogen_half_saturation" => Dict("a" => 1, "b" => 1),
        ),
        "linear_mortality" => 8e-7 / second,
        "alpha" => 0.1953 / day,
    ),
    zoo_args=Dict(
        "volumes" => Dict(
            "min_volume" => 10, "max_volume" => 100, "splitting" => "linear_splitting"
        ),
        "allometry" => Dict("maximum_predation_rate" => Dict("a" => 1, "b" => 1)),
        "linear_mortality" => 8e-7 / second,
        "holling_half_saturation" => 5.0,
        "quadratic_mortality" => 1e-6 / second,
    ),
    palatability_args=Dict(
        "P" => Dict("optimum_predator_prey_ratio" => 0, "protection" => 0),
        "Z" => Dict("optimum_predator_prey_ratio" => 10, "protection" => 1),
    ),
    assimilation_efficiency_args=Dict(
        "P" => Dict("can_be_eaten" => 1, "can_eat" => 0, "assimilation_efficiency" => 0),
        "Z" => Dict("can_be_eaten" => 0, "can_eat" => 1, "assimilation_efficiency" => 0.32),
    ),
    bgc_args=Dict(
        "detritus_remineralization" => 0.1213 / day, "mortality_export_fraction" => 0.5
    ),
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
)
    phyto_args["n"] = n_phyto
    zoo_args["n"] = n_zoo

    # compute emergent parameters
    defined_parameters = Dict("P" => phyto_args, "Z" => zoo_args)

    if isnothing(palatability_matrix)
        defined_parameters["P"]["palatability"] = palatability_args["P"]
        defined_parameters["Z"]["palatability"] = palatability_args["Z"]
    end

    if isnothing(assimilation_efficiency_matrix)
        defined_parameters["P"]["assimilation_efficiency"] = assimilation_efficiency_args["P"]
        defined_parameters["Z"]["assimilation_efficiency"] = assimilation_efficiency_args["Z"]
    end

    emergent_parameters = compute_allometric_parameters(defined_parameters)

    # combine emergent parameters with remaining user defined parameters
    parameters = NamedTuple(
        Symbol(k) => v for (k, v) in merge(bgc_args, emergent_parameters)
    )

    if !isnothing(palatability_matrix)
        parameters["palatability_matrix"] = palatability_matrix
    end

    if !isnothing(assimilation_efficiency_matrix)
        parameters["assimilation_efficiency_matrix"] = assimilation_efficiency_matrix
    end

    # create tracer functions
    plankton_array = vcat(
        [Symbol("P$i") for i in 1:n_phyto], [Symbol("Z$i") for i in 1:n_zoo]
    )
    tracers = Dict(
        "N" => nutrient_dynamics(plankton_array), "D" => detritus_dynamics(plankton_array)
    )
    for i in 1:n_phyto
        name = "P$i"
        tracers[name] = phyto_dynamics(plankton_array, name)
    end
    for i in 1:n_zoo
        name = "Z$i"
        tracers[name] = zoo_dynamics(plankton_array, name)
    end

    # return Oceananigans.Biogeochemistry object
    return define_tracer_functions(parameters, tracers)
end

end # module
