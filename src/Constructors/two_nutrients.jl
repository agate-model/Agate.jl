module Constructors

include("Biogeochemistry.jl")
include("Parameters.jl")
include("Tracers.jl")

using .Biogeochemistry
using .Parameters
using .Tracers

using Oceananigans.Units

export construct_size_structured_P_Z_POC_DOC_DIN_PO4

"""
Construct an instance of xyz model.

Note that if non-default `*_dynamics` expressions are passed, the relevant `*_args` also
need to be specified.

# Arguments
- `n_phyto`: number of phytoplankton to include in the model
- `n_zoo`: number of zooplankton to include in the model
- `DIC_dynamics`: 
- `PO4_dynamics`:
- `DIN_dynamics`:
- `POC_dynamics`:
- `DOC_dynamics`:
- `phyto_dynamics`: expression describing how phytoplankton grow
- `zoo_dynamics`: expression describing how zooplankton grow
- `phyto_args`: Dictionary of phytoplankton parameters
- `zoo_args`: Dictionary of zooplankton parameters
- `interaction_args`: Dictionary of arguments from which a palatability and assimilation
   efficiency matrix between all plankton can be computed
- `bgc_args`: biogeochemistry parameters related to nutrient and detritus
- `palatability_matrix`: optional palatability matrix passed as a NamedArray, if provided
   then `palatability_args` are ignored
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as a
   NamedArray, if provided then `assimilation_args` are ignored
"""
function construct_size_structured_P_Z_POC_DOC_DIN_PO4(;
    n_phyto=2,
    n_zoo=2,
    DIC_dynamics = DIC_typical,
    PO4_dynamics=PO4_fixed_quota,
    DIN_dynamics=DIN_fixed_quota,
    POC_dynamics=POC_typical,
    DOC_dynamics=DOC_typical,
    phyto_dynamics=phytoplankton_growth_fixed_quota_nitrogen_and_phosphorus,
    zoo_dynamics=zooplankton_growth_simplified,
    phyto_args=Dict(
        "diameters" =>
            Dict("min_diameter" => 2, "max_diameter" => 10, "splitting" => "log_splitting"),
        "allometry" => Dict(
            "maximum_growth_rate" => Dict("a" => 2 / day, "b" => -0.15),
            "half_saturation_DIN" => Dict("a" => 0.17, "b" => 0.27),
            "half_saturation_PO4" => Dict("a" => 0.17, "b" => 0.27),
        ),
        "linear_mortality" => 8e-7 / second,
        "alpha" => 0.1953 / day,
    ),
    zoo_args=Dict(
        "diameters" => Dict(
            "min_diameter" => 20,
            "max_diameter" => 100,
            "splitting" => "linear_splitting",
        ),
        "allometry" =>
            Dict("maximum_predation_rate" => Dict("a" => 30.84 / day, "b" => -0.16)),
        "linear_mortality" => 8e-7 / second,
        "holling_half_saturation" => 5.0,
        "quadratic_mortality" => 1e-6 / second,
    ),
    interaction_args=Dict(
        "P" => Dict(
            "can_eat" => 0,
            "can_be_eaten" => 1,
            "optimum_predator_prey_ratio" => 0,
            "protection" => 0,
            "specificity" => 0,
            "assimilation_efficiency" => 0,
        ),
        "Z" => Dict(
            "can_eat" => 1,
            "can_be_eaten" => 0,
            "optimum_predator_prey_ratio" => 10,
            "protection" => 1,
            "specificity" => 0.3,
            "assimilation_efficiency" => 0.32,
        ),
    ),
    bgc_args=Dict(
        "POC_remineralization" => 0.1213 / day, 
        "DOC_remineralization" => 0.1213 / day, 
        "mortality_export_fraction" => 0.5,
        "nitrogen_to_carbon" => 1,
        "phosphorus_to_carbon" => 1,
    ),
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
)
    phyto_args["n"] = n_phyto
    zoo_args["n"] = n_zoo

    # compute emergent parameters
    defined_parameters = Dict("P" => phyto_args, "Z" => zoo_args)

    if isnothing(palatability_matrix)
        defined_parameters["P"]["palatability"] = Dict(
            k => interaction_args["P"][k] for
            k in ["can_eat", "optimum_predator_prey_ratio", "protection", "specificity"]
        )
        defined_parameters["Z"]["palatability"] = Dict(
            k => interaction_args["Z"][k] for
            k in ["can_eat", "optimum_predator_prey_ratio", "protection", "specificity"]
        )
    end

    if isnothing(assimilation_efficiency_matrix)
        defined_parameters["P"]["assimilation_efficiency"] = Dict(
            k => interaction_args["P"][k] for
            k in ["can_eat", "can_be_eaten", "assimilation_efficiency"]
        )
        defined_parameters["Z"]["assimilation_efficiency"] = Dict(
            k => interaction_args["Z"][k] for
            k in ["can_eat", "can_be_eaten", "assimilation_efficiency"]
        )
    end

    emergent_parameters = compute_allometric_parameters(defined_parameters)

    if !isnothing(palatability_matrix)
        if !(size(palatability_matrix) == (n_phyto + n_zoo, n_phyto + n_zoo))
            throw(
                ArgumentError(
                    "palatability_matrix must have size $((n_phyto+n_zoo, n_phyto+n_zoo))"
                ),
            )
        end
        emergent_parameters["palatability_matrix"] = palatability_matrix
    end

    if !isnothing(assimilation_efficiency_matrix)
        if !(size(assimilation_efficiency_matrix) == (n_phyto + n_zoo, n_phyto + n_zoo))
            throw(
                ArgumentError(
                    "assimilation_efficiency_matrix must have size $((n_phyto+n_zoo, n_phyto+n_zoo))",
                ),
            )
        end

        emergent_parameters["assimilation_efficiency_matrix"] =
            assimilation_efficiency_matrix
    end

    # combine emergent parameters with remaining user defined parameters
    parameters = NamedTuple(
        Symbol(k) => v for (k, v) in merge(bgc_args, emergent_parameters)
    )

    # create tracer functions
    plankton_array = vcat(
        [Symbol("P$i") for i in 1:n_phyto], [Symbol("Z$i") for i in 1:n_zoo]
    )
    tracers = Dict(
        "DIC" => DIC_dynamics(plankton_array), 
        "DIN" => DIN_dynamics(plankton_array), 
        "PO4" => PO4_dynamics(plankton_array), 
        "POC" => POC_dynamics(plankton_array),
        "DOC" => DOC_dynamics(plankton_array),
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
