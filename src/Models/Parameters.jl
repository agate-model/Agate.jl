module Parameters

using NamedArrays

export compute_darwin_parameters

# TODO: the DARWIN emergent functions should eventually be defined here - use dummy for now
include("../../examples/emergent_4P/emergent_functions.jl")
emergent_max_growth_rate_f = dummy_emergent_growth
emergent_max_predation_rate_f = dummy_emergent_predation_rate
emergent_nitrogen_half_saturation_f = dummy_emergent_nitrogen_half_saturation
emergent_palatability_f = dummy_emergent_palat
emergent_assimilation_efficiency_f = dummy_emergent_assimilation_efficiency

# parameters the user has to pass to compute_darwin_parameters
EXPECTED_EMERGENT_PARAMS = [
    "growth_a",
    "growth_b",
    "nitrogen_half_saturation_a",
    "nitrogen_half_saturation_b",
    "predation_rate_a",
    "predation_rate_b",
    "optimum_predator_prey_ratio",
    "protection",
]
EXPECTED_VOLUME_PARAMS = ["n", "min_volume", "max_volume", "splitting"]

"""
Generate `n` volumes for each `plankton` species (e.g., "P", "Z" or "cocco") and for each
species-volume combination compute:
    - maximum growth rate
    - maximum predation rate
    - nitrogen half saturation
    - assimilation efficiency
    - palatability

The `plankton` Dictionary keys have to contain argument names of the emergent functions as
well as `["n", "min_volume", "max_volume", "splitting"]` and should look something like:
    ```
    plankton = Dict(
        <species name> => Dict(
            "n" => <value>,
            "min_volume" => <value>,
            "max_volume" => <value>,
            "splitting" => <one of "linear_splitting" or "log_splitting">,
            <function arg name> => <value>,
            ...
        ),
        <species name> => Dict(...),
    )
    ```
"""
function compute_darwin_parameters(plankton::Dict)

    # validate that `plankton` has the required parameters for each plankton species
    for (plankton_name, params) in plankton
        for p in vcat(EXPECTED_EMERGENT_PARAMS, EXPECTED_VOLUME_PARAMS)
            if !(p ∈ keys(params))
                throw(ArgumentError("$plankton_name parameter dictionary missing $p"))
            end
        end
    end

    # generate n volumes for the species and compute emergent parameters
    parameters_with_volume = split_size_parameters(plankton)

    emergent_results = Dict()
    for (f, name) in zip(
        [
            emergent_max_growth_rate_f,
            emergent_max_predation_rate_f,
            emergent_nitrogen_half_saturation_f,
        ],
        ["maximum_growth_rate", "maximum_predation_rate", "nitrogen_half_saturation"],
    )
        emergent_results[name] = emergent_1D_array(parameters_with_volume, f)
    end

    for (f, name) in zip(
        [emergent_palatability_f, emergent_assimilation_efficiency_f],
        [["maximum_growth_rate", "maximum_predation_rate", "nitrogen_half_saturation"]],
    )
        emergent_results[name] = emergent_2D_array(parameters_with_volume, f)
    end

    # add remaining parameters (except those used in emergent computations)
    # we are assuming each species has the same parameters - so just get the first name
    species = collect(keys(parameters_with_volume))
    for param_name in keys(parameters_with_volume[species[1]])
        if !(param_name ∈ EXPECTED_EMERGENT_PARAMS)
            # Collect values across species in an array
            values = [
                parameters_with_volume[species][param_name] for
                species in keys(parameters_with_volume)
            ]
            species_names = collect(keys(parameters_with_volume))

            # Store as a NamedArray
            emergent_results[param_name] = NamedArray(
                values, (species_names,), ("Species",)
            )
        end
    end

    return emergent_results
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

"""
Generate `n` volumes for each `plankton` species (e.g., "P", "Z" or "cocco") and return
Dictionary of parameters for each species-volume combination.

The returned Dictionary keys contain `volume` as well as any other species specific parameters
passed in the `plankton` Dictionary.
"""
function split_size_parameters(plankton::Dict)
    parameters_with_volume = Dict()

    for (plankton_name, params) in plankton
        n = params["n"]
        # expect the splitting function is defined within this module
        splitting_function = getfield(Parameters, Symbol(params["splitting"]))
        volumes = splitting_function(n, params["min_volume"], params["max_volume"])

        for p in EXPECTED_VOLUME_PARAMS
            delete!(params, p)
        end

        for i in 1:n
            parameters_with_volume["$plankton_name$i"] = merge(
                Dict("volume" => volumes[i]), params
            )
        end
    end

    return parameters_with_volume
end

"""
Calculate value for each `plankton` species using `func`, returning a 1D NamedArray.
The functions passed here are expected to return one of:
    - maximum growth rate
    - maximum predation rate
    - nitrogen half saturation

The `plankton` Dictionary keys have to contain argument names of the function to apply to it
and should be of the form `Dict(<species name> => Dict(<arg name> => <val>, ...), ... )`
"""
function emergent_1D_array(plankton::Dict, func::Function)
    # Get species names
    species_names = collect(keys(plankton))

    # first argname is self --> skip
    argnames = Base.method_argnames.(methods(func))[1][2:end]

    # Calculate emergent values for each species
    emergent_values = [
        func([plankton[name][String(key)] for key in argnames]...) for name in species_names
    ]

    # Create a NamedArray with species names as row labels
    emergent_array = NamedArray(emergent_values, (species_names,), ("Species",))

    return emergent_array
end

"""
Calculate value for every pair of `plankton` species using `func`, returning a 2D NamedArray.
The functions passed here are expected to return one of:
    - assimilation efficiency
    - palatability

The `plankton` Dictionary keys have to contain argument names of the function to apply to it
and should be of the form `Dict(<species name> => Dict(<arg name> => <val>, ...), ... )`
"""
function emergent_2D_array(plankton::Dict, func::Function)
    # Extract predator and prey names
    predator_names = collect(keys(plankton))
    prey_names = collect(keys(plankton))

    # Initialize a NamedArray to hold values
    values = zeros(Float64, length(predator_names), length(prey_names))
    names = (predator=predator_names, prey=prey_names)

    # Populate the NamedArray with calculated values
    for (i, pred_name) in enumerate(predator_names)
        for (j, prey_name) in enumerate(prey_names)
            prey_data = plankton[prey_name]
            predator_data = plankton[pred_name]

            # Pass prey and predator data dictionaries to the function
            values[i, j] = func(prey_data, predator_data)
        end
    end

    # Create and return a NamedArray with the values and names
    return NamedArray(values, names)
end

end #module
