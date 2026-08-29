using Agate
using Test
using ForwardDiff

using Agate.Library.Allometry:
    consumer_assimilation_matrix_axes, palatability_matrix_allometric_axes,
    resolve_diameter_indexed_vector
using Agate.Library.Nutrients:
    frank_tnorm, liebig_minimum, normalized_droop_limitation, quota_uptake_regulation
using Agate.Library.Photosynthesis: geider_light_response, smith_light_limitation
using Agate.Library.Predation: holling_type_ii, preferential_predation_loss

@testset "Library" begin
    @test holling_type_ii(1.0, 1.0) == 0.5
    @test preferential_predation_loss(0.2, 1.0, 0.5, 0.1, 0.2, 0.8) ≈ 1 / 150
end

@testset "Allometry accepts realized diameter tuples" begin
    diameters = (1.0, 2.0)
    @test resolve_diameter_indexed_vector(Float64, diameters, (2,), 3.0; default=0.0) == [0.0, 3.0]
    @test palatability_matrix_allometric_axes(
        Float64, diameters;
        optimum_predator_prey_ratio=[0.0, 2.0],
        specificity=[0.0, 1.0],
        protection=[0.0, 0.0],
        consumer_indices=(2,),
        prey_indices=(1, 2),
    ) == [1.0 0.5]
    @test consumer_assimilation_matrix_axes(
        Float64; assimilation_efficiency=[0.2, 0.8], consumer_indices=(2,), prey_indices=(1, 2)
    ) == [0.8 0.8]
end

@testset "Library scalar genericity" begin
    T = Float32

    @test Agate.Library.Nutrients.monod_limitation(T(1), T(0.5)) isa T
    @test frank_tnorm(T(0.2), T(0.4)) isa T
    @test frank_tnorm(T(0.2), T(0.4); sharpness=50.0) isa T
    @test smith_light_limitation(T(50), T(0.1), T(1)) isa T
    @test Agate.Library.Mortality.linear_loss(T(1), T(0.1)) isa T
    @test Agate.Library.Predation.holling_type_ii(T(1), T(0.5)) isa T
    @test Agate.Library.Remineralization.linear_remineralization(T(1), T(0.1)) isa T
    @test Agate.Library.Temperature.q10_temperature_factor(T(10), T(2)) isa T
    @test Agate.Library.Temperature.q10_temperature_factor(T(30), T(2), T(20)) == T(2)
    geider_response = geider_light_response(T(100), T(2e-6), T(2e-5), T(0.02))
    @test geider_response isa T
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

    gradient = ForwardDiff.gradient(x -> frank_tnorm(x[1], x[2]), [0.5, 0.5])
    @test all(isfinite, gradient)
    @test gradient[1] ≈ gradient[2]
end

@testset "Droop quota kernels" begin
    qmin, qmax = 0.1, 0.2
    reference = 1.0

    @test normalized_droop_limitation(0.05, reference, qmin, qmax) == 0.0
    @test normalized_droop_limitation(qmin, reference, qmin, qmax) == 0.0
    @test normalized_droop_limitation(qmax, reference, qmin, qmax) == 1.0
    @test normalized_droop_limitation(0.3, reference, qmin, qmax) == 1.0
    @test normalized_droop_limitation(0.15, reference, qmin, qmax) ≈ 2 / 3
    @test normalized_droop_limitation(0.0, 0.0, qmin, qmax) == 0.0

    @test quota_uptake_regulation(0.05, reference, qmin, qmax, 2.0) == 1.0
    @test quota_uptake_regulation(qmin, reference, qmin, qmax, 2.0) == 1.0
    @test quota_uptake_regulation(qmax, reference, qmin, qmax, 2.0) == 0.0
    @test quota_uptake_regulation(0.3, reference, qmin, qmax, 2.0) == 0.0
    @test quota_uptake_regulation(0.15, reference, qmin, qmax, 2.0) ≈ 0.75
    @test quota_uptake_regulation(0.0, 0.0, qmin, qmax, 2.0) == 0.0

    @test 0.0 <= normalized_droop_limitation(0.15, reference, qmin, qmax) <= 1.0
    @test 0.0 <= quota_uptake_regulation(0.15, reference, qmin, qmax, 2.0) <= 1.0
    derivative = ForwardDiff.derivative(
        internal -> normalized_droop_limitation(internal, reference, qmin, qmax), 0.15
    )
    @test 0 < derivative < Inf

    T = Float32
    @test normalized_droop_limitation(T(0.15), one(T), T(qmin), T(qmax)) isa T
    @test quota_uptake_regulation(T(0.15), one(T), T(qmin), T(qmax), T(2)) isa T
end
