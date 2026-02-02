using Agate
using Test

using Agate.Library.Predation:
    HollingTypeII,
    IdealizedPredationLoss

@testset "Library" begin
    @test HollingTypeII(1.0)(1.0) == 0.5

    loss = IdealizedPredationLoss(0.1, 0.2)(1.0, 0.5)
    @test loss > 0
end
