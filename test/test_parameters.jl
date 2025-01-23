using Agate

using Oceananigans.Units

include(joinpath("..", "src", "Library", "Library.jl"))

using .Library.Allometry
using .Library.Predation

@testset "Models.Parameters" begin
    @testset "create_allometric_parameters" begin

        # include the 0 parameters to make it comparable to the N2P2ZD example
        defined_parameters = Dict(
            "P" => Dict(
                "n" => 2,
                "diameters" => Dict(
                    "min_diameter" => 2,
                    "max_diameter" => 10,
                    "splitting" => "log_splitting",
                ),
                "allometry" => Dict(
                    "maximum_growth_rate" => Dict("a" => 2 / day, "b" => -0.15),
                    "nutrient_half_saturation" => Dict("a" => 0.17, "b" => 0.27),
                    "maximum_predation_rate" => Dict("a" => 0, "b" => 0),
                ),
                "palatability" => Dict(
                    "can_eat" => 0,
                    "optimum_predator_prey_ratio" => 0,
                    "protection" => 0,
                    "specificity" => 0,
                ),
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
                "diameters" => Dict(
                    "min_diameter" => 20,
                    "max_diameter" => 100,
                    "splitting" => "linear_splitting",
                ),
                "allometry" => Dict(
                    "maximum_growth_rate" => Dict("a" => 0, "b" => 0),
                    "nutrient_half_saturation" => Dict("a" => 0, "b" => 0),
                    "maximum_predation_rate" => Dict("a" => 30.84 / day, "b" => -0.16),
                ),
                "palatability" => Dict(
                    "can_eat" => 1,
                    "optimum_predator_prey_ratio" => 10,
                    "protection" => 1,
                    "specificity" => 0.3,
                ),
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

        emergent_parameters = compute_allometric_parameters(defined_parameters)

        # compare against hand computed `parameters` in examples
        include(joinpath("..", "examples", "N2P2ZD", "tracers.jl"))

        plankton_order = ["P1", "P2", "Z1", "Z2"]

        for (key, emerge_params) in emergent_parameters
            # start with arrays of values
            if !(
                key in
                ["assimilation_efficiency_matrix", "palatability_matrix", "diameters"]
            )
                @test all(
                    isapprox.(
                        parameters[Symbol(key)],
                        [emerge_params[p] for p in plankton_order],
                        rtol=0.01,
                    ),
                ) || println(
                    "Test failed for parameter: $(key). Expected: $(parameters[Symbol(key)]), Got: $([emerge_params[p] for p in plankton_order])",
                )
                # matrices of values -> compare row at a time
            elseif !(key == "diameters")
                for (i, p) in enumerate(plankton_order)
                    emergent_row = emerge_params[p, :]
                    true_row = parameters[Symbol(replace(key, "_matrix" => ""))][i, :]
                    @test all(
                        isapprox.(
                            true_row, [emergent_row[p] for p in plankton_order], rtol=0.01
                        ),
                    ) || println(
                        "Test failed for parameter: $(key), row: $(i). Expected: $(true_row), Got: $([emergent_row[p] for p in plankton_order])",
                    )
                end
            end
        end
    end
end
