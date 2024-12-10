
@testset "Models.Tracers" begin
    # N2P2ZD model
    # TODO: use functions in Models library
    include(joinpath("..", "examples", "emergent_4P", "tracers.jl"))
    model = N2P2ZD()
    P1 = 0.01
    P2 = 0.01
    Z1 = 0.05
    Z2 = 0.05
    N = 7.0
    D = 1
    PAR = 100

    @test !iszero(model(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
    @test !iszero(model(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D, PAR))
end
