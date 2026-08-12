using Agate
using Test
using ForwardDiff

using Agate.Library.Nutrients: liebig_minimum, smooth_liebig_minimum
using Agate.Library.Photosynthesis: liebig_nutrient_limitation, smooth_liebig_nutrient_limitation
using Agate.Library.Predation: holling_type_ii, idealized_predation_loss

@testset "Library" begin
    @test holling_type_ii(1.0, 1.0) == 0.5

    loss = idealized_predation_loss(1.0, 0.5, 0.1, 0.2)
    @test loss > 0
end

@testset "Library scalar genericity" begin
    T = Float32

    @test Agate.Library.Nutrients.monod_limitation(T(1), T(0.5)) isa T
    @test Agate.Library.Nutrients.smooth_liebig_minimum(T(0.2), T(0.4)) isa T
    @test Agate.Library.Photosynthesis.smooth_liebig_nutrient_limitation(
        (T(1), T(2)), (T(0.5), T(1)), one(T)
    ) isa T
    @test Agate.Library.Photosynthesis.smith_light_limitation(T(50), T(0.1), T(1)) isa T
    @test Agate.Library.Mortality.linear_loss(T(1), T(0.1)) isa T
    @test Agate.Library.Predation.holling_type_ii(T(1), T(0.5)) isa T
    @test Agate.Library.Remineralization.linear_remineralization(T(1), T(0.1)) isa T
    @test Agate.Library.Temperature.q10_temperature_factor(T(10), T(2)) isa T
end
