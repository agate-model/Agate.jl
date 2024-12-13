using Agate

@testset "Models.Constructor" begin
    # N2P2ZD model
    include(joinpath("..", "examples", "emergent_4P", "tracers.jl"))

    # TODO: eventually use constructor here
    include(joinpath("..", "examples", "N2P2ZD", "tracers.jl"))

    model = N2P2ZD()
    # model_emergent = N2P2ZD_emergent()
    N2P2ZD_emergent = construct_NPZD_instance()
    model_emergent = N2P2ZD_emergent()

    P1 = 0.01
    P2 = 0.01
    Z1 = 0.05
    Z2 = 0.05
    N = 7.0
    D = 1
    PAR = 100

    @test !iszero(model_emergent(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model_emergent(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model_emergent(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model_emergent(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model_emergent(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model_emergent(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))

    @test model_emergent(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
        model(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
    @test model_emergent(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
        model(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
    @test model_emergent(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
        model(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
    @test model_emergent(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
        model(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
    @test model_emergent(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
        model(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
    @test model_emergent(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR) ==
        model(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR)
end
