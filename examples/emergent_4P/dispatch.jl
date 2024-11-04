include("emergent_functions.jl")
include("helper_functions.jl")

# Sample plankton data with required keys
plankton = Dict(
    "species1" => Dict(
        "volume" => 1.0,
        "volume_a" => 1.0,
        "volume_b" => 1.0,
        "pred" => 2.0,
        "optimum_predator_prey_ratio" => 10,
        "protection" => 0,
    ),
    "species2" => Dict(
        "volume" => 1.5,
        "volume_a" => 1.0,
        "volume_b" => 1.0,
        "pred" => 1.8,
        "optimum_predator_prey_ratio" => 10,
        "protection" => 0,
    ),
)

"""

"""
function emergent_analysis(plankton, func, params)
    result = nothing
    try
        result = emergent_1D_array(plankton, func, params)
        if ndims(result_1D) == 1
            println("1D")
        else
            error("Unexpected output type")
        end
    catch
    end
    try
        result = emergent_2D_array(plankton, func, params)
        if ndims(result_2D) == 2
            println("2D")
        else
            error("Unexpected output type")
        end
    catch
    end

    return result
end

emergent_analysis(
    plankton, dummy_emergent_palat, ["volume", "optimum_predator_prey_ratio", "protection"]
)
emergent_analysis(
    plankton, dummy_emergent_predation_rate, ["volume_a", "volume_b", "volume"]
)
