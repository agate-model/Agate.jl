"""
Module to compute size-dependent plankton parameters.
"""

module Parameters

using DataStructures: DefaultDict

using Agate.Library.Allometry
using Agate.Library.Predation

export compute_allometric_parameters

emergent_palatability_f = allometric_palatability_unimodal_protection
emergent_assimilation_efficiency_f = assimilation_efficiency_emergent_binary

"""
    compute_allometric_parameters(plankton::Dict) -> Tuple(Dict, Array)

This function:
    - generates `n` diameters for each `plankton` group using either a linear or a log
      splitting scale
    - optionally computes emergent parameters (allometric functions, assimilation matrix,
      palatability matrix) for each group-diameter combination
    - reshapes any other group specific parameters to Array of length `n` (repeating the
      parameter value n times to match it in length to the emergent parameters)

All parameters are returned as: `Dict(<parameter> => <Array of values>, ....)` along with
an Array of the plankton order in which the parameter values were created and processed
(e.g., `["Z1", ...., "Z<n>", "P1", ..., "P<n>"]`).

# Arguments
- `plankton`: a Dictionary of plankton groups' specific parameters of the form:
       `Dict(<group name> => Dict(<parameter> => <value>, ....), ...)`

    The Dictionary for each group (e.g., "P", "Z" or "cocco") has to contain at least the
    keys "n" and "diameters", which are either an array of values or a dictionary of the
    following form:
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

    The Dictionary for each group can optionally include the following keys and values:
        - key: "allometry" (compute diameter dependent parameters)
            - value: Dictionary of the form
                `Dict(<param name> => Dict("a" => <value>, "b" => <value>), ...)`
        - key: "palatability" (generate a palatability matrix)
            - value: Dictionary with "can_eat", "optimum_predator_prey_ratio", "protection",
              "specificity" keys
        - key: "assimilation_efficiency" (generate assimilation efficiency matrix)
            - value: Dictionary with "can_eat", "can_be_eaten", "assimilation_efficiency" keys
    For example:
        ```
        Dict(
            "P" => Dict(
                ...,
                "allometry" => Dict(
                    "maximum_growth_rate" => Dict("a" => 2.3148e-5, "b" => -0.15),
                    "nutrient_half_saturation" => Dict("a" => 0.17, "b" => 0.27),
                ),
                ...
            ),
            "Z" => Dict(
                ...,
                "palatability" => Dict(
                    "can_eat" => 1,
                    "optimum_predator_prey_ratio" => 10,
                    "protection" => 1,
                    "specificity" => 0.3,
                ),
                "assimilation_efficiency" => Dict(
                    "can_be_eaten" => 0,
                    "can_eat" => 1,
                    "assimilation_efficiency" => 0.32
                ),
                ...
            )
        )
        ```
"""
function compute_allometric_parameters(plankton::Dict)

    # intermediate representations for parameters that are output as matrices
    intermediate_palatability = DefaultDict{AbstractString,Array}([])
    intermediate_assimilation = DefaultDict{AbstractString,Array}([])
    # final outputs here
    results = DefaultDict{AbstractString,Array}([])

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

        # 2. compute allometric functions (if any specified by user)
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
                    # NOTE: expect here that in all other cases `value` is a single number
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
Apply function `func` to every every pair of `plankton` groups, returning a 2D Array. The
expected use is for calculating predator-prey dynamics (palatability and assimilation
efficiency matrices).

The `plankton` Dictionary keys have to contain argument names of the function to calculate
and be of the form `Dict(<function arg name> => <Array of values>, ...)`

# Example
```julia
function add_ab(prey, pred)
    return prey["a"] + pred["b]
end

plankton = Dict(
    "a" => [1, 2, 3],
    "b" => [9, 8, 7],
)

values_matrix = emergent_2D_array(plankton, add_ab)
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

end #module
