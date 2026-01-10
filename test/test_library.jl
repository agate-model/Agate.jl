using Agate
using Test

using Agate.Library.Predation:
    AssimilationPreyParameters,
    AssimilationPredatorParameters,
    holling_type_2,
    predation_loss_idealized,
    assimilation_efficiency_emergent_binary

@testset "Library" begin
    @test holling_type_2(1.0, 1.0) == 0.5

    loss = predation_loss_idealized(0.1, 0.2, 1.0, 0.5)
    @test loss > 0

    prey = AssimilationPreyParameters(true)
    predator = AssimilationPredatorParameters{Float32}(true, 0.32f0)
    @test assimilation_efficiency_emergent_binary(prey, predator) == 0.32f0

    prey2 = AssimilationPreyParameters(false)
    @test assimilation_efficiency_emergent_binary(prey2, predator) == 0.0f0
end
