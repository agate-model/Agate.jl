using Agate
using Test

using Agate.Library.Predation:
    holling_type_ii,
    idealized_predation_loss

@testset "Library" begin
    @test holling_type_ii(1.0, 1.0) == 0.5

    loss = idealized_predation_loss(1.0, 0.5, 0.1, 0.2)
    @test loss > 0
end
