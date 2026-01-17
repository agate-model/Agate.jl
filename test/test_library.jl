using Agate
using Test

using Agate.Library.Predation:
    AssimilationPreyParameters,
    AssimilationPredatorParameters,
    HollingTypeII,
    IdealizedPredationLoss,
    EmergentAssimilationEfficiencyBinary

@testset "Library" begin
    @test HollingTypeII(1.0)(1.0) == 0.5

    loss = IdealizedPredationLoss(0.1, 0.2)(1.0, 0.5)
    @test loss > 0

    prey = AssimilationPreyParameters(true)
    predator = AssimilationPredatorParameters{Float32}(true, 0.32f0)
    @test EmergentAssimilationEfficiencyBinary()(prey, predator) == 0.32f0

    prey2 = AssimilationPreyParameters(false)
    @test EmergentAssimilationEfficiencyBinary()(prey2, predator) == 0.0f0
end
