using Agate
using Test
using ForwardDiff

using Agate.Library.Nutrients: liebig_minimum, smooth_liebig_minimum
using Agate.Library.Photosynthesis: liebig_nutrient_limitation, smooth_liebig_nutrient_limitation
using Agate.Library.Predation: holling_type_ii, idealized_predation_loss
using Agate.Tendencies: TendencyConfig

@testset "Library" begin
    @test holling_type_ii(1.0, 1.0) == 0.5

    loss = idealized_predation_loss(1.0, 0.5, 0.1, 0.2)
    @test loss > 0
end
