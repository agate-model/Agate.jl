using Agate
using Agate.Models.Parameters: split_size_parameters, emergent_1D_array, emergent_2D_array

include(joinpath("..", "examples", "emergent_4P", "emergent_functions.jl"))

@testset "Models.Parameters" begin
    @testset "size splitting" begin
        defined_parameters = Dict(
            "P" => Dict(
                "n" => 2,
                "min_volume" => 1,
                "max_volume" => 10,
                "splitting" => "log_splitting",
                "growth_a" => 1,
                "assimilation_efficiency" => 0,
            ),
            "Z" => Dict(
                "n" => 2,
                "min_volume" => 1,
                "max_volume" => 100,
                "splitting" => "linear_splitting",
                "growth_a" => 0,
                "assimilation_efficiency" => 1,
            ),
        )

        intermediate_parameters = split_size_parameters(defined_parameters)
        @test isapprox(intermediate_parameters["P1"]["volume"], 1)
        @test isapprox(intermediate_parameters["Z2"]["volume"], 100)
        @test isapprox(intermediate_parameters["Z2"]["growth_a"], 0)
        @test isapprox(intermediate_parameters["Z2"]["assimilation_efficiency"], 1)
    end

    @testset "1D NamedArray" begin
        plankton = Dict(
            "P1" => Dict("growth_a" => 1, "growth_b" => 0, "volume" => 1),
            "Z1" => Dict("growth_a" => 0, "growth_b" => 0, "volume" => 10.0),
        )

        result = emergent_1D_array(plankton, dummy_emergent_growth)
        @test isapprox(result["P1"], 7.190e-6)
        @test isapprox(result["Z1"], 0)
    end

    @testset "2D NamedArray" begin
        plankton = Dict(
            "P1" =>
                Dict("volume" => 1, "optimum_predator_prey_ratio" => 0, "protection" => 0),
            "P2" =>
                Dict("volume" => 10, "optimum_predator_prey_ratio" => 0, "protection" => 0),
            "Z1" => Dict(
                "volume" => 10, "optimum_predator_prey_ratio" => 10, "protection" => 1
            ),
            "Z2" => Dict(
                "volume" => 100, "optimum_predator_prey_ratio" => 10, "protection" => 1
            ),
        )

        result = emergent_2D_array(plankton, dummy_emergent_palat)
        @test isapprox(result["P1", "P2"], 0)
        @test isapprox(result["Z1", "P2"], 0.3)
        @test isapprox(result["Z2", "P2"], 1)
    end

end
