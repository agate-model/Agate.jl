"""
Module to compute and format size-dependent plankton parameters.
"""

module Parameters

using DataStructures: DefaultDict

using Agate.Library.Allometry
using Agate.Library.Predation

export compute_allometric_parameters, create_size_structured_params

emergent_palatability_f = allometric_palatability_unimodal_protection
emergent_assimilation_efficiency_f = assimilation_efficiency_emergent_binary

"""
    compute_allometric_parameters(plankton::Dict) -> Tuple(Dict, Array)

Compute allometric (diameter-dependant) plankton parameters.

# Arguments
- `plankton`: a Dictionary of plankton groups' specific parameters of the form:
       `Dict(<group name> => Dict(<parameter> => <value>, ....), ...)`
    The Dictionary for each group (e.g., "P", "Z" or "cocco") has to contain at least the
    keys "n" and "diameters", specifying number of plankton diameters to generate and how
    (min/max values and whether to use linear or log splitting), e.g.:
        ```
        Dict(
            "P" => Dict(
                "n" => 2,
                "diameters" =>
                    Dict("min_diameter" => 1, "max_diameter" => 10, "splitting" => "log_splitting"),
                ...
                ),
            "Z" => ...
        )
        ```

!!! info

    The Dictionary for each plankton group can also include the following keys and values:
        - key: "allometry" (compute diameter dependent parameters)
            - value: Dictionary of the form
                `Dict(<param name> => Dict("a" => <value>, "b" => <value>), ...)`
        - key: "palatability" (generate a palatability matrix)
            - value: Dictionary with "can_eat", "optimum_predator_prey_ratio", "protection",
              "specificity" keys
        - key: "assimilation_efficiency" (generate assimilation efficiency matrix)
            - value: Dictionary with "can_eat", "can_be_eaten", "assimilation_efficiency" keys
        - key: "*" (any other Strings)
             - value: Float specifying group level parameters (e.g., mortality rates)

!!! info

    This function:
    - generates `n` diameters for each `plankton` group
    - computes emergent parameters (allometric functions, assimilation matrix,
      palatability matrix) for each group-diameter combination, if these are specified as keys
      in the input Dictionary
    - reshapes any other group specific parameters to Array of length `n` (repeating the
      parameter value n times to match it in length to the emergent parameters)

    All parameters are returned as: `Dict(<parameter> => <Array of values>, ....)` along with
    an Array of the plankton order in which the parameter values were created and processed
    (e.g., `["Z1", ...., "Z<n>", "P1", ..., "P<n>"]`).

!!! tip

    See test.test_parameters `create_allometric_parameters` testset for example input.
"""
function compute_allometric_parameters(plankton::Dict)

    # intermediate representations for parameters that are output as matrices
    intermediate_palatability = DefaultDict{AbstractString,Array}([])
    intermediate_assimilation = DefaultDict{AbstractString,Array}([])
    # final outputs here
    results = DefaultDict{AbstractString,Array{Real}}([])

    plankton_names = []
    for (plankton_name, params) in plankton
        n = params["n"]
        plankton_names = vcat(plankton_names, ["$plankton_name$i" for i in 1:n])

        # 1. compute diameters (unless already specified)
        if isa(params["diameters"], Dict)
            # get the appropriate splitting function from the Parameters module
            splitting_function = getfield(
                Parameters, Symbol(params["diameters"]["splitting"])
            )
            diameters = splitting_function(
                params["diameters"]["min_diameter"], params["diameters"]["max_diameter"], n
            )
        else
            if length(params["diameters"]) != n
                throw(ArgumentError("diameters array must have length $(n)"))
            end
            diameters = params["diameters"]
        end
        results["diameters"] = vcat(results["diameters"], diameters)

        # 2. compute allometric functions (if "allometry" specified by user in params)
        if "allometry" ∈ keys(params)
            for (param, args) in params["allometry"]
                values = [
                    allometric_scaling_power(args["a"], args["b"], d) for d in diameters
                ]
                results[param] = vcat(results[param], values)
            end
        end

        # 3. reshape remaining parameters
        for (param, value) in params
            if !(param ∈ ["n", "diameters", "allometry"])

                # NOTE: for palatability and assimilation_efficiency, `value` is a Dictionary
                # store the values in intermediate form to use later when creating matrices
                if param ∈ ["palatability", "assimilation_efficiency"]
                    for (arg, val) in value
                        values = repeat([val], n)
                        if param == "palatability"
                            intermediate_palatability[arg] = vcat(
                                intermediate_palatability[arg], values
                            )
                        elseif param == "assimilation_efficiency"
                            intermediate_assimilation[arg] = vcat(
                                intermediate_assimilation[arg], values
                            )
                        end
                    end
                else
                    # NOTE: expect here that in all other cases `value` is a Float specifying
                    # group level parameters
                    results[param] = vcat(results[param], repeat([value], n))
                end
            end
        end
    end

    # 4. palatability & assimilation efficiency matrices (if specified by user)
    if !(isempty(intermediate_palatability))
        intermediate_palatability["diameters"] = results["diameters"]
        results["palatability_matrix"] = emergent_2D_array(
            intermediate_palatability, emergent_palatability_f
        )
    end

    if !(isempty(intermediate_assimilation))
        intermediate_assimilation["diameters"] = results["diameters"]
        results["assimilation_efficiency_matrix"] = emergent_2D_array(
            intermediate_assimilation, emergent_assimilation_efficiency_f
        )
    end

    return results, plankton_names
end

"""
Log splitting function to generate a set of diameters.
"""
function log_splitting(min_diameter::Real, max_diameter::Real, n::Int)
    log_min = log10(min_diameter)
    log_max = log10(max_diameter)
    log_step = (log_max - log_min) / (n - 1)
    return [10^(log_min + i * log_step) for i in 0:(n - 1)]
end

"""
Linear splitting function to generate a set of diameters.
"""
function linear_splitting(min_diameter::Real, max_diameter::Real, n::Int)
    linear_step = (max_diameter - min_diameter) / (n - 1)
    return [min_diameter + i * linear_step for i in 0:(n - 1)]
end

"""
    emergent_2D_array(plankton, func)

Apply function `func` to every every pair of `plankton` groups, returning a 2D Array. The
expected use is for calculating predator-prey dynamics (palatability and assimilation
efficiency matrices).

# Arguments
- `plankton`: Dictionary of values to apply `func` to, of the form
   `Dict(<function arg name> => <Array of values>, ...)`
- `func`: Function

# Example
```julia
function add_ab(prey, pred)
    return prey["a"] + pred["b]
end

plankton_vals = Dict(
    "a" => [1, 2, 3],
    "b" => [9, 8, 7],
)

values_matrix = emergent_2D_array(plankton_vals, add_ab)
```
"""
function emergent_2D_array(plankton, func)
    # n is the number of plankton in the system - populate NxN matrix
    n = length(plankton["diameters"])
    plankton_matrix = zeros(Float64, n, n)

    # each plankton can be a predator or prey, retrieve their respective args
    arg_names = keys(plankton)
    for pred_name in 1:n
        predator_data = Dict(arg => plankton[arg][pred_name] for arg in arg_names)
        for prey_name in 1:n
            prey_data = Dict(arg => plankton[arg][prey_name] for arg in arg_names)
            # pass prey and predator data dictionaries to the function
            plankton_matrix[pred_name, prey_name] = func(prey_data, predator_data)
        end
    end

    return plankton_matrix
end

"""
    create_size_structured_params(;
        n_plankton=Dict("P" => 2, "Z" => 2),
        diameters=Dict(
            "P" =>
                Dict("min_diameter" => 2, "max_diameter" => 10, "splitting" => "log_splitting"),
            "Z" => Dict(
                "min_diameter" => 20,
                "max_diameter" => 100,
                "splitting" => "linear_splitting",
            ),
        ),
        plankton_args=nothing,
        interaction_args=nothing,
        bgc_args=nothing,
        palatability_matrix=nothing,
        assimilation_efficiency_matrix=nothing,
    )

Create a NamedTuple of parameters from which an Oceananigans.Biogeochemistry model can be instantiated.

# Arguments
- `n_plankton`: Dictionary of the number of plankton to include in the model by group
- `diameters`: Dictionary which specifies for each plankton group how to compute diameters
   or gives a list of values to use
- `plankton_args`: Dictionary of plankton parameters for each group
- `interaction_args`: Dictionary of arguments from which a palatability and assimilation
   efficiency matrix between all plankton can be computed
- `bgc_args`: Dictionary of constant parameters used in growth functions (i.e., not size
    dependant plankton parameters as well as biogeochemistry parameters related to nutrient
    and detritus
- `palatability_matrix`: optional palatability matrix passed as an Array, if provided
    then `interaction_args` are not used to compute this
- `assimilation_efficiency_matrix`: optional assimilation efficiency matrix passed as an
    Array, if provided then `interaction_args` are not used to compute this

!!! warning

    Wherever a parameter is defined for only some plankton groups, its value is set to 0 for
    the other plankton groups.

!!! tip

    See Agate.Models.NiPiZD.Constructor for example plankton_args, interaction_args and bgc_args.

"""
function create_size_structured_params(;
    n_plankton=Dict("P" => 2, "Z" => 2),
    diameters=Dict(
        "P" =>
            Dict("min_diameter" => 2, "max_diameter" => 10, "splitting" => "log_splitting"),
        "Z" => Dict(
            "min_diameter" => 20,
            "max_diameter" => 100,
            "splitting" => "linear_splitting",
        ),
    ),
    plankton_args=nothing,
    interaction_args=nothing,
    bgc_args=nothing,
    palatability_matrix=nothing,
    assimilation_efficiency_matrix=nothing,
)
    # ====================================================================================
    # 1. Ensure that all plankton groups have the same parameters
    #   - gather all parameters across groups
    #   - when param not specified for a group, add it and set value to 0
    # ====================================================================================

    # NOTE: Julia passes dictionaries by reference not value
    # use deepcopy here to avoid mutating the passed dictionaries outside this fucntion
    plankton_args_copy = deepcopy(plankton_args)

    # first handle allometry
    # - get all allometric parameters across groups (e.g., Ps and Zs)
    all_allometric_params = union(
        vcat([collect(keys(args["allometry"])) for args in values(plankton_args_copy)]...)
    )
    # - loop across groups, if allometric param not defined, set it to 0
    for name in keys(plankton_args_copy)
        for param in all_allometric_params
            if !(param ∈ keys(plankton_args_copy[name]["allometry"]))
                plankton_args_copy[name]["allometry"][param] = Dict("a" => 0.0, "b" => 0.0)
            end
        end
    end

    # second, handle non allometric params (same as above, exclude allometry)
    all_other_params = setdiff(
        union(vcat([collect(keys(args)) for args in values(plankton_args_copy)]...)),
        ["allometry"],
    )
    for name in keys(plankton_args_copy)
        for param in all_other_params
            if !(param ∈ keys(plankton_args_copy[name]))
                plankton_args_copy[name][param] = 0.0
            end
        end
    end

    # ====================================================================================
    # 2. Create dictionary that matches `compute_allometric_parameters` expected input
    # ====================================================================================

    # add remaining paramaters that need to compute allometric values
    for name in keys(plankton_args_copy)
        plankton_args_copy[name]["n"] = n_plankton[name]
        plankton_args_copy[name]["diameters"] = diameters[name]
    end

    defined_parameters = plankton_args_copy

    # ====================================================================================
    # 3. Handle matrix inputs (these are either provided by the user or will be computed)
    # ====================================================================================

    if isnothing(palatability_matrix)
        for name in keys(plankton_args_copy)
            defined_parameters[name]["palatability"] = Dict(
                k => interaction_args[name][k] for
                k in ["can_eat", "optimum_predator_prey_ratio", "protection", "specificity"]
            )
        end
    end

    if isnothing(assimilation_efficiency_matrix)
        for name in keys(plankton_args_copy)
            defined_parameters[name]["assimilation_efficiency"] = Dict(
                k => interaction_args[name][k] for
                k in ["can_eat", "can_be_eaten", "assimilation_efficiency"]
            )
        end
    end

    # ====================================================================================
    # 4. Compute allometric parameters
    # ====================================================================================

    # NOTE: function also reshapes non-emergent plankton parameters
    emergent_parameters, plankton_names = compute_allometric_parameters(defined_parameters)

    # ====================================================================================
    # 5. Clean up
    # ====================================================================================

    # if user provided matrix inputs, add them here (check if correct size)
    n = sum(values(n_plankton))
    if !isnothing(palatability_matrix)
        if !(size(palatability_matrix) == (n, n))
            throw(ArgumentError("palatability_matrix must have size $((n, n))"))
        end
        emergent_parameters["palatability_matrix"] = palatability_matrix
    end

    if !isnothing(assimilation_efficiency_matrix)
        if !(size(assimilation_efficiency_matrix) == (n, n))
            throw(ArgumentError("assimilation_efficiency_matrix must have size $((n, n))"))
        end

        emergent_parameters["assimilation_efficiency_matrix"] =
            assimilation_efficiency_matrix
    end

    # append information that need access to at model object instantiation
    for (name, n) in n_plankton
        bgc_args["n_$name"] = n
    end

    # combine emergent parameters with remaining bgc parameters
    parameters = NamedTuple(
        Symbol(k) => v for (k, v) in merge(bgc_args, emergent_parameters)
    )

    return parameters, plankton_names
end

end #module
