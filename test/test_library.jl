using Agate
using Test
using ForwardDiff

using Agate.Library.Nutrients: FrankMinimum, frank_minimum, liebig_minimum
using Agate.Library.Photosynthesis: frank_nutrient_limitation, liebig_nutrient_limitation
using Agate.Library.Predation: holling_type_ii, idealized_predation_loss
using Agate.Tendencies: TendencyConfig, nutrient_coupling, phytoplankton_tendency

@testset "Library" begin
    @test holling_type_ii(1.0, 1.0) == 0.5

    loss = idealized_predation_loss(1.0, 0.5, 0.1, 0.2)
    @test loss > 0
end

@testset "Library scalar genericity" begin
    T = Float32

    @test Agate.Library.Nutrients.monod_limitation(T(1), T(0.5)) isa T
    @test frank_minimum(T(0.2), T(0.4)) isa T
    @test frank_minimum(T(0.2), T(0.4); sharpness=50.0) isa T
    @test frank_nutrient_limitation((T(1), T(2)), (T(0.5), T(1)), one(T)) isa T
    @test frank_nutrient_limitation(
        (T(1), T(2)), (T(0.5), T(1)), one(T); sharpness=25.0
    ) isa T
    @test Agate.Library.Photosynthesis.smith_light_limitation(T(50), T(0.1), T(1)) isa T
    @test Agate.Library.Mortality.linear_loss(T(1), T(0.1)) isa T
    @test Agate.Library.Predation.holling_type_ii(T(1), T(0.5)) isa T
    @test Agate.Library.Remineralization.linear_remineralization(T(1), T(0.1)) isa T
    @test Agate.Library.Temperature.q10_temperature_factor(T(10), T(2)) isa T
end

@testset "Frank minimum" begin
    @test frank_minimum(0.5, 1.0) ≈ 0.5
    @test frank_minimum(0.0, 0.5) ≈ 0.0
    @test frank_minimum(1.0, 1.0) ≈ 1.0
    @test frank_minimum(0.2, 0.8) ≈ frank_minimum(0.8, 0.2)
    @test frank_minimum((0.2, 1.0, 0.8)) ≈ frank_minimum(0.2, 0.8)
    @test isnan(frank_minimum(NaN, 0.5))

    low_sharpness = frank_minimum(0.5, 0.5; sharpness=25)
    high_sharpness = frank_minimum(0.5, 0.5; sharpness=100)
    @test low_sharpness < high_sharpness < liebig_minimum(0.5, 0.5)

    s = 5.0
    q = exp(-s)
    expected = log(1 + ((q^0.3 - 1) * (q^0.7 - 1)) / (q - 1)) / log(q)
    @test frank_minimum(0.3, 0.7; sharpness=s) ≈ expected

    @test frank_nutrient_limitation((1.0,), (1.0,), 1.0) ≈ 0.5

    gradient = ForwardDiff.gradient(x -> FrankMinimum()(x[1], x[2]), [0.5, 0.5])
    @test all(isfinite, gradient)
    @test gradient[1] ≈ gradient[2]

    liebig_config = TendencyConfig(;
        growth=:smith, organic_cycling=:simple_detritus, nutrient_limitation=:liebig
    )
    @test liebig_config.nutrient_limitation isa Agate.Library.Nutrients.LiebigMinimum

    frank_config = TendencyConfig(;
        growth=:smith,
        organic_cycling=:simple_detritus,
        nutrient_limitation=FrankMinimum(25),
    )
    @test frank_config.nutrient_limitation isa FrankMinimum
    @test frank_config.nutrient_limitation.sharpness == 25
end

@testset "Frank tendency plumbing" begin
    config = TendencyConfig(;
        growth=:geider,
        organic_cycling=:dom_pom,
        nutrient_limitation=FrankMinimum(50),
        nutrients=(
            nutrient_coupling(
                :DIN,
                :half_saturation_DIN;
                stoichiometry=:nitrogen_to_carbon,
                remineralization=(
                    (:DON, :DON_remineralization), (:PON, :PON_remineralization)
                ),
            ),
            nutrient_coupling(
                :PO4,
                :half_saturation_PO4;
                stoichiometry=:phosphorus_to_carbon,
                remineralization=(
                    (:DOP, :DOP_remineralization), (:POP, :POP_remineralization)
                ),
            ),
        ),
    )

    bgc = Agate.Models.DARWIN.construct()
    frank_tendency = phytoplankton_tendency(config; plankton_idx=3)

    DIN = bgc.parameters.half_saturation_DIN[3]
    PO4 = bgc.parameters.half_saturation_PO4[3]
    args = (
        10.0,
        DIN,
        PO4,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.01,
        0.01,
        100.0,
    )

    frank = frank_tendency(bgc, 0, 0, 0, 0, args...)
    liebig = bgc(Val(:P_1), 0, 0, 0, 0, args...)

    @test isfinite(frank)
    @test frank < liebig
end
