using Agate
using NamedArrays

@testset "Models.Constructor" begin
    @testset "N2P2ZD model" begin
        # N2P2ZD model defined using low level syntax
        include(joinpath("..", "examples", "N2P2ZD", "tracers.jl"))
        model = N2P2ZD()

        # N2P2ZD model constructed from emergent parameters
        N2P2ZD_constructed = construct_size_structured_NPZD()
        model_constructed = N2P2ZD_constructed()

        P1 = 0.01
        P2 = 0.01
        Z1 = 0.05
        Z2 = 0.05
        N = 7.0
        D = 1
        PAR = 100

        @test !iszero(model_constructed(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
        @test !iszero(model_constructed(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))

        @test model_constructed(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
            model(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
        @test model_constructed(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
            model(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
        @test model_constructed(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
            model(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
        @test model_constructed(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
            model(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
        @test model_constructed(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
            model(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
        @test model_constructed(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
            model(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
    end

    @testset "User defined matrices" begin
        names = ["P", "Z"]
        wrong_size_matrix = NamedArray(zeros(Float64, 2, 2), (predator=names, prey=names))
        @test_throws ArgumentError construct_size_structured_NPZD(
            palatability_matrix=wrong_size_matrix
        )
        @test_throws ArgumentError construct_size_structured_NPZD(
            assimilation_efficiency_matrix=wrong_size_matrix
        )

        # doesn't throw error if dimensions are correct
        names = ["P1", "P2", "Z1", "Z2"]
        correct_size_matrix = NamedArray(zeros(Float64, 4, 4), (predator=names, prey=names))
        new_model = construct_size_structured_NPZD(;
            palatability_matrix=correct_size_matrix,
            assimilation_efficiency_matrix=correct_size_matrix,
        )
    end
end
