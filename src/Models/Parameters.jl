module Parameters

using DataStructures: DefaultDict
using NamedArrays

export compute_darwin_parameters

# TODO: the DARWIN emergent functions should eventually be defined here - use dummy for now
include(joinpath("..", "..", "examples", "emergent_4P", "emergent_functions.jl"))
emergent_palatability_f = dummy_emergent_palat
emergent_assimilation_efficiency_f = dummy_emergent_assimilation_efficiency

"""
THIS FUNCTION IS A PLACEHOLDER -- TO BE UPDATED!
"""
function allometry_f(param, a, b, volume)
    if param == "maximum_growth_rate"
        return dummy_emergent_growth(a, b, volume)
    elseif param == "nitrogen_half_saturation"
        return dummy_emergent_nitrogen_half_saturation(a, b, volume)
    elseif param == "maximum_predation_rate"
        return dummy_emergent_predation_rate(a, b, volume)
    end
end

"""
    compute_darwin_parameters(plankton::Dict) -> Dict

Generate `n` volumes for each `plankton` species (e.g., "P", "Z" or "cocco") and for each
species-volume combination compute:
    - maximum_growth_rate
    - maximum_predation_rate
    - nitrogen_half_saturation
    - assimilation_efficiency
    - palatability
All computed and other species specific parameters are returned as a Dictionary of the form:
`Dict(<parameter> => <NamedArray of values>, ....)`.

# Arguments
- `plankton`: a Dictionary of plankton species specific parameters of the form:
       `Dict(<species name> => Dict(<parameter> => <value>, ....), ...)`
   the species all have to have the same parameters and these have to include at least the
   parameters in `EXPECTED_EMERGENT_PARAMS` and `EXPECTED_VOLUME_PARAMS`
"""
function compute_darwin_parameters(plankton::Dict)

    # intermediate representations
    intermediate_palatability = DefaultDict{AbstractString,NamedArray}(NamedArray([], []))
    intermediate_assimilation = DefaultDict{AbstractString,NamedArray}(NamedArray([], []))
    # final outputs here
    results = DefaultDict{AbstractString,NamedArray}(NamedArray([], []))

    for (plankton_name, params) in plankton
        n = params["n"]

        plankton_names = ["$plankton_name$i" for i in 1:n]

        # 1. compute plankton volumes
        splitting_function = getfield(Parameters, Symbol(params["volumes"]["splitting"]))
        volumes = splitting_function(
            params["volumes"]["min_volume"], params["volumes"]["max_volume"], n
        )
        results["volumes"] = vcat(results["volumes"], NamedArray(volumes, plankton_names))

        # 2. compute allometric functions if specified
        if "allometry" ∈ keys(results)
            for (param, args) in params["allometry"]
                # TODO: once `allometry_f` is updated, it should not have `param` as argument
                values = [allometry_f(param, args["a"], args["b"], v) for v in volumes]
                results[param] = vcat(result[name], NamedArray(values, plankton_names))
            end
        end

        # 3. reshape remaining parameters
        # NOTE: for palatability and assimilation_efficiency, `value` is a Dictionary
        for (param, value) in params
            if !(param ∈ ["n", "volumes", "allometry"])
                if param == "palatability"
                    for (arg, val) in value
                        intermediate_palatability[arg] = vcat(
                            intermediate_palatability[arg],
                            NamedArray(repeat([val], n), plankton_names),
                        )
                    end
                elseif param == "assimilation_efficiency"
                    for (arg, val) in value
                        intermediate_assimilation[arg] = vcat(
                            intermediate_assimilation[arg],
                            NamedArray(repeat([val], n), plankton_names),
                        )
                    end
                else
                    # in this case val is a single value
                    results[param] = vcat(
                        results[param], NamedArray(repeat([value], n), plankton_names)
                    )
                end
            end
        end
    end

    # 4. TODO: palatability & assimilation efficiency
    if !(isempty(intermediate_palatability))
        # println(intermediate_palatability)
        intermediate_palatability["volumes"] = results["volumes"]
        results["palatability_matrix"] = emergent_2D_array(
            intermediate_palatability, emergent_palatability_f
        )
    end

    if !(isempty(intermediate_assimilation))
        # println(intermediate_assimilation)
        intermediate_assimilation["volumes"] = results["volumes"]
        results["assimilation_efficiency_matrix"] = emergent_2D_array(
            intermediate_assimilation, emergent_assimilation_efficiency_f
        )
    end

    return results
end

"""
Log splitting function to generate a set of volumes.
"""
function log_splitting(min_volume::Real, max_volume::Real, n::Int)
    log_min = log10(min_volume)
    log_max = log10(max_volume)
    log_step = (log_max - log_min) / (n - 1)
    return [10^(log_min + i * log_step) for i in 0:(n - 1)]
end

"""
Linear splitting function to generate a set of volumes.
"""
function linear_splitting(min_volume::Real, max_volume::Real, n::Int)
    linear_step = (max_volume - min_volume) / (n - 1)
    return [min_volume + i * linear_step for i in 0:(n - 1)]
end

# TODO: update docstring
"""
Apply function `func` to every every pair of `plankton` groups, returning a 2D NamedArray.

The `plankton` Dictionary keys have to contain argument names of the function to calculate
and be of the form `Dict(<function arg name> => NamedArray(<values>, <names>), ...)`

The function takes in a pair of NamedArray that capture the required values...

# Example
```julia
function check(prey, pred) return prey["a"] + pred["b] end
plankton = Dict(
    "a" => NamedArray([1, 2, 3], ["A1", "B2", "C3"]),
    "b" => NamedArray([9, 8, 7], ["A1", "B2", "C3"]),
)
values_matrix = emergent_2D_array(plankton, check)
```
"""
function emergent_2D_array(plankton, func)
    plankton_names = names(plankton["volumes"])[1]
    arg_names = keys(plankton)
    values = zeros(Float64, length(plankton_names), length(plankton_names))
    plankton_matrix = NamedArray(values, (predator=plankton_names, prey=plankton_names))

    # Q: is there a better way to populate this ?!
    # Populate the NamedArray with calculated values
    for pred_name in plankton_names
        for prey_name in plankton_names
            prey_data = Dict(arg => plankton[arg][prey_name] for arg in arg_names)
            predator_data = Dict(arg => plankton[arg][pred_name] for arg in arg_names)
            # Pass prey and predator data dictionaries to the function
            plankton_matrix[pred_name, prey_name] = func(prey_data, predator_data)
        end
    end

    # Create and return a NamedArray with the values and names
    return plankton_matrix
end

end #module
