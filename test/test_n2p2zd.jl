using Agate.Models.Dynamic: expression_check
using Oceananigans.Biogeochemistry:
    AbstractContinuousFormBiogeochemistry,
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields

@testset "N2P2ZD" begin
    using Oceananigans.Units
    # N2P2ZD model
    include("../examples/N2P2ZD/model.jl")
    model = N2P2ZD()
    P1 = 0.01
    P2 = 0.01
    Z1 = 0.05
    Z2 = 0.05
    N = 7.0
    D = 1
    model(Val(:N), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D)
    model(Val(:D), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D)
    model(Val(:P1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D)
    model(Val(:P2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D)
    model(Val(:Z1), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D)
    model(Val(:Z2), 0, 0, 0, 0, P1, P2, Z1, Z2, N, D)
end
