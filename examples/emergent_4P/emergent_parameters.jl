using NamedArrays

include("emergent_functions.jl")
include("helper_functions.jl")

# Example plankton data
defined_parameters = Dict(
    "P1" => Dict(
        "volume" => 1,
        "growth_a" => 1,
        "growth_b" => 1,
        "protection" => 0,
        "opt_ratio" => 0,
    ),
    "P2" => Dict(
        "volume" => 10,
        "growth_a" => 1,
        "growth_b" => 1,
        "protection" => 0,
        "opt_ratio" => 0,
    ),
    "Z1" => Dict(
        "volume" => 10,
        "growth_a" => 0,
        "growth_b" => 0,
        "protection" => 1,
        "opt_ratio" => 10,
    ),
    "Z2" => Dict(
        "volume" => 100,
        "growth_a" => 0,
        "growth_b" => 0,
        "protection" => 1,
        "opt_ratio" => 10,
    ),
)

emergent_functions = Dict(
    "growth_rate" => :(dummy_emergent_growth(growth_a, growth_b, volume)),
    "palatability" => :(dummy_emergent_palat(
        prey_volume, predator_volume, optimum_predator_prey_ratio, protection
    )),
    "predation_rate" => :(dummy_emergent_predation_rate(volume_a, volume_b, volume)),
    "nitrogen_half_saturation" => :(dummy_emergent_nitrogen_half_saturation(volume_a, volume_b, volume)),
)

#then pass to multiple dispatch