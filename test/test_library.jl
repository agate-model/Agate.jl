using Agate
using Test
using ForwardDiff

using Agate.Library.Nutrients: frank_tnorm, liebig_minimum
using Agate.Library.Photosynthesis:
    frank_nutrient_limitation, geider_light_limitation, geider_light_response,
    liebig_nutrient_limitation
using Agate.Library.Predation: holling_type_ii, idealized_predation_loss, preferential_predation_loss

@testset "Library" begin
    @test holling_type_ii(1.0, 1.0) == 0.5
    @test idealized_predation_loss(2.0, 0.5, 0.1, 1.0) ≈ 0.04

    @test preferential_predation_loss(1.0, 0.5, 0.1, 0.2, 0.8) ≈ 1 / 30
end

@testset "Library scalar genericity" begin
    T = Float32

    @test Agate.Library.Nutrients.monod_limitation(T(1), T(0.5)) isa T
    @test frank_tnorm(T(0.2), T(0.4)) isa T
    @test frank_tnorm(T(0.2), T(0.4); sharpness=50.0) isa T
    @test frank_nutrient_limitation((T(1), T(2)), (T(0.5), T(1)), one(T)) isa T
    @test frank_nutrient_limitation(
        (T(1), T(2)), (T(0.5), T(1)), one(T); sharpness=25.0
    ) isa T
    @test Agate.Library.Photosynthesis.smith_light_limitation(T(50), T(0.1), T(1)) isa T
    @test Agate.Library.Mortality.linear_loss(T(1), T(0.1)) isa T
    @test Agate.Library.Predation.holling_type_ii(T(1), T(0.5)) isa T
    @test Agate.Library.Remineralization.linear_remineralization(T(1), T(0.1)) isa T
    @test Agate.Library.Temperature.q10_temperature_factor(T(10), T(2)) isa T
    @test Agate.Library.Temperature.q10_temperature_factor(T(30), T(2), T(20)) == T(2)
    geider_response = geider_light_response(T(100), T(2e-6), T(2e-5), T(0.02))
    @test geider_response isa T
    @test geider_light_limitation(T(100), T(2e-6), T(2e-5), T(0.02)) == T(2e-5) * geider_response
end

@testset "Frank t-norm" begin
    @test frank_tnorm(0.5, 1.0) ≈ 0.5
    @test frank_tnorm(0.0, 0.5) ≈ 0.0
    @test frank_tnorm(1.0, 1.0) ≈ 1.0
    @test frank_tnorm(0.2, 0.8) ≈ frank_tnorm(0.8, 0.2)
    @test frank_tnorm((0.2, 1.0, 0.8)) ≈ frank_tnorm(0.2, 0.8)
    @test isnan(frank_tnorm(NaN, 0.5))

    low_sharpness = frank_tnorm(0.5, 0.5; sharpness=25)
    high_sharpness = frank_tnorm(0.5, 0.5; sharpness=100)
    @test low_sharpness < high_sharpness < liebig_minimum(0.5, 0.5)

    s = 5.0
    q = exp(-s)
    expected = log(1 + ((q^0.3 - 1) * (q^0.7 - 1)) / (q - 1)) / log(q)
    @test frank_tnorm(0.3, 0.7; sharpness=s) ≈ expected

    @test frank_nutrient_limitation((1.0,), (1.0,), 1.0) ≈ 0.5

    gradient = ForwardDiff.gradient(x -> frank_tnorm(x[1], x[2]), [0.5, 0.5])
    @test all(isfinite, gradient)
    @test gradient[1] ≈ gradient[2]

    resources = (0.5, 0.5)
    half_saturations = (0.5, 0.5)
    frank = frank_nutrient_limitation(
        resources, half_saturations, 1.0; sharpness=50
    )
    liebig = liebig_nutrient_limitation(resources, half_saturations, 1.0)
    @test isfinite(frank)
    @test frank < liebig

end
