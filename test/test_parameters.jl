using Agate

using Oceananigans.Units

# TODO: this will eventually be replaced with the actual emergent functions
include(joinpath("..", "examples", "emergent_4P", "emergent_functions.jl"))

@testset "Models.Parameters" begin
    @testset "create_darwin_parameters" begin

        # include the 0 parameters to make it comparable to the N2P2ZD example
        defined_parameters = Dict(
            "P" => Dict(
                "n" => 2,
                "volumes" => Dict(
                    "min_volume" => 1,
                    "max_volume" => 10,
                    "splitting" => "log_splitting",
                ),
                "allometry" => Dict(
                    "maximum_growth_rate" => Dict("a" => 1, "b" => 1),
                    "nitrogen_half_saturation" => Dict("a" => 1, "b" => 1),
                    "maximum_predation_rate" => Dict("a" => 0, "b" => 0),
                ),
                "palatability" =>
                    Dict("optimum_predator_prey_ratio" => 0, "protection" => 0),
                "assimilation_efficiency" => Dict(
                    "can_be_eaten" => 1,
                    "can_eat" => 0,
                    "assimilation_efficiency" => 0,
                ),
                "linear_mortality" => 8e-7 / second,
                "holling_half_saturation" => 0,
                "quadratic_mortality" => 0,
                "alpha" => 0.1953 / day,
            ),
            "Z" => Dict(
                "n" => 2,
                "volumes" => Dict(
                    "min_volume" => 10,
                    "max_volume" => 100,
                    "splitting" => "linear_splitting",
                ),
                "allometry" => Dict(
                    "maximum_growth_rate" => Dict("a" => 0, "b" => 0),
                    "nitrogen_half_saturation" => Dict("a" => 0, "b" => 0),
                    "maximum_predation_rate" => Dict("a" => 1, "b" => 1),
                ),
                "palatability" =>
                    Dict("optimum_predator_prey_ratio" => 10, "protection" => 1),
                "assimilation_efficiency" => Dict(
                    "can_be_eaten" => 0,
                    "can_eat" => 1,
                    "assimilation_efficiency" => 0.32,
                ),
                "linear_mortality" => 8e-7 / second,
                "holling_half_saturation" => 5.0,
                "quadratic_mortality" => 1e-6 / second,
                "alpha" => 1e-99 / day,
            ),
        )

        emergent_parameters = compute_darwin_parameters(defined_parameters)

        # compare against hand computer `parameters` in examples
        include(joinpath("..", "examples", "N2P2ZD", "tracers.jl"))

        plankton_order = ["P1", "P2", "Z1", "Z2"]

        for (key, emerge_params) in emergent_parameters
            # start with arrays of values
            if !(key ∈ ["assimilation_efficiency_matrix", "palatability_matrix", "volumes"])
                @test all(
                    parameters[Symbol(key)] .== [emerge_params[p] for p in plankton_order]
                )
                # matrices of values -> compare row at a time
            elseif !(key == "volumes")
                for (i, p) in enumerate(plankton_order)
                    emergent_row = emerge_params[p, :]
                    true_row = parameters[Symbol(replace(key, "_matrix" => ""))][i, :]
                    @test all(true_row .== [emergent_row[p] for p in plankton_order])
                end
            end
        end
    end
end
