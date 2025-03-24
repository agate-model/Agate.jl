module Utils

using Agate.Models.Parameters

export create_params_dict

"""
Create a dictionary of parameters to pass to `Agate.Models.Biogeochemistry.define_tracer_functions`.

Used for models with phytoplankton and zooplankton (NiPiZD, DARWIN).

Wherever a parameter is defined for only one plankton group, its value is set to 0 for the
other group. This way all the returned parameter Arrays are of same length (n_phyto + n_zoo).

# Arguments
- `n_phyto`: number of phytoplankton to include in the model
- `n_zoo`: number of zooplankton to include in the model
- `phyto_diameters`: dictionary from which `n_phyto` diameters can be computed or a list of
    values to use
- `zoo_diameters`: dictionary from which `zoo` diameters can be computed or a list of
    values to use
- `phyto_args`: Dictionary of phytoplankton parameters, for default values see
    `Agate.Models.Constructors.DEFAULT_PHYTO_ARGS`
- `zoo_args`: Dictionary of zooplankton parameters, for default values see
    `Agate.Models.Constructors.DEFAULT_ZOO_ARGS`
- `interaction_args`: Dictionary of arguments from which a palatability and assimilation
   efficiency matrix between all plankton can be computed, for default values see
    `Agate.Models.Constructors.DEFAULT_INTERACTION_ARGS`
- `bgc_args`: Dictionary of constant parameters used in growth functions (i.e., not size
    dependant plankton parameters as well as biogeochemistry parameters related to nutrient
    and detritus, for default values see `Agate.Models.Constructors.DEFAULT_CONSTANT_ARGS`
- `palatability_matrix`: optional palatability matrix passed as an Array, if provided
    then `interaction_args` are not used to compute this
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as an
    Array, if provided then `interaction_args` are not used to compute this
"""
function create_params_dict(;
    n_phyto=2,
    n_zoo=2,
    phyto_diameters=Dict(
        "min_diameter" => 2, "max_diameter" => 10, "splitting" => "log_splitting"
    ),
    zoo_diameters=Dict(
        "min_diameter" => 20, "max_diameter" => 100, "splitting" => "linear_splitting"
    ),
    phyto_args=nothing,
    zoo_args=nothing,
    interaction_args=nothing,
    bgc_args=nothing,
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
)
    # ====================================================================================
    # 1. Make sure that phyto_args and zoo_args have all the same parameters to ensure the
    #    output parameter vectors all have the same size (n_phyto + n_zoo)
    #    - set value of parameter to 0 in those cases
    # ====================================================================================

    # NOTE: Julia passes dictionaries by reference not value
    # use deepcopy here to avoid mutating the passed dictionaries outside this fucntion
    phyto_args_copy = deepcopy(phyto_args)
    zoo_args_copy = deepcopy(zoo_args)

    # first handle allometry
    all_allometric_params = union(
        keys(phyto_args["allometry"]), keys(zoo_args["allometry"])
    )
    for param in all_allometric_params
        if !(param ∈ keys(phyto_args["allometry"]))
            phyto_args_copy["allometry"][param] = Dict("a" => 0, "b" => 0)
        end
        if !(param ∈ keys(zoo_args["allometry"]))
            zoo_args_copy["allometry"][param] = Dict("a" => 0, "b" => 0)
        end
    end

    # then handle non allometric params
    all_other_params = setdiff(union(keys(phyto_args), keys(zoo_args)), ["allometry"])
    for param in all_other_params
        if !(param ∈ keys(phyto_args))
            phyto_args_copy[param] = 0
        end
        if !(param ∈ keys(zoo_args))
            zoo_args_copy[param] = 0
        end
    end

    # ====================================================================================
    # 2. Create dictionary that matches `compute_allometric_parameters` expected input
    # ====================================================================================

    # add remaining paramaters that need to compute allometric values
    phyto_args_copy["n"] = n_phyto
    phyto_args_copy["diameters"] = phyto_diameters
    zoo_args_copy["n"] = n_zoo
    zoo_args_copy["diameters"] = zoo_diameters

    defined_parameters = Dict("P" => phyto_args_copy, "Z" => zoo_args_copy)

    # ====================================================================================
    # 3. Handle matrix inputs (these are either provided by the user or will be computed)
    # ====================================================================================

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

    # ====================================================================================
    # 4. Compute allometric parameters
    #   - if user provided matrix inputs, add them here (check if correct size)
    # ====================================================================================

    # NOTE: function also reshapes non-emergent plankton parameters
    emergent_parameters, plankton_names = compute_allometric_parameters(defined_parameters)

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

    # ====================================================================================
    # 5. Clean up
    # ====================================================================================

    # append information that need access to at instantiation
    bgc_args["n_phyto"] = n_phyto
    bgc_args["n_zoo"] = n_zoo

    # combine emergent parameters with remaining bgc parameters
    parameters = NamedTuple(
        Symbol(k) => v for (k, v) in merge(bgc_args, emergent_parameters)
    )

    return parameters, plankton_names
end

end # module
