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

@testset "AD-friendly nutrient limitation" begin
    values = (0.2, 0.5, 0.8)
    @test liebig_minimum(values) == 0.2
    @test smooth_liebig_minimum(values; sharpness=100.0) ≈ 0.2 atol=1e-3

    equal_values = (0.5, 0.5)
    hard_slope = ForwardDiff.derivative(x -> liebig_minimum((x, 0.5)), 0.5)
    smooth_slope = ForwardDiff.derivative(x -> smooth_liebig_minimum((x, 0.5)), 0.5)
    @test isfinite(hard_slope)
    @test smooth_slope ≈ 0.5 atol=1e-8

    resources = (1.0, 0.25)
    half_saturations = (1.0, 0.25)
    hard = liebig_nutrient_limitation(resources, half_saturations, 1.0)
    smooth = smooth_liebig_nutrient_limitation(resources, half_saturations, 1.0)
    @test smooth <= hard
    @test isfinite(smooth)

    config = TendencyConfig(; growth=:smith,
                              organic_cycling=:simple_detritus,
                              nutrient_limitation=:smooth_liebig)
    @test config isa TendencyConfig
end
