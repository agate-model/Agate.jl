using Agate
using Agate.Models: NiPiZD
using Agate.Models.Parameters: PhytoPFTParameters
using Agate.Models.NiPiZD.Tracers: nutrient_geider_light, phytoplankton_geider_light

using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, biogeochemical_drift_velocity

@testset "Models.NiPiZD constructors" begin
    bgc_type = NiPiZD.construct()
    model = bgc_type()

    @test required_biogeochemical_tracers(model) == (:N, :D, :Z1, :Z2, :P1, :P2)

    P1 = 0.01
    P2 = 0.01
    Z1 = 0.05
    Z2 = 0.05
    N = 7.0
    D = 1.0
    PAR = 100.0

    tracer_vals(sym) =
        if sym === :P1
            P1
        elseif sym === :P2
            P2
        elseif sym === :Z1
            Z1
        elseif sym === :Z2
            Z2
        elseif sym === :N
            N
        else
            D
        end
    ordered = [tracer_vals(s) for s in required_biogeochemical_tracers(model)]

    @test isfinite(model(Val(:N), 0, 0, 0, 0, ordered..., PAR))
    @test isfinite(model(Val(:D), 0, 0, 0, 0, ordered..., PAR))
    @test isfinite(model(Val(:P1), 0, 0, 0, 0, ordered..., PAR))
    @test isfinite(model(Val(:Z1), 0, 0, 0, 0, ordered..., PAR))

    @testset "User-defined interaction matrices" begin
        wrong = zeros(Float64, 2, 2)
        @test_throws ArgumentError NiPiZD.instantiate(bgc_type; palatability_matrix=wrong)
        @test_throws ArgumentError NiPiZD.instantiate(
            bgc_type; assimilation_efficiency_matrix=wrong
        )

        correct = zeros(Float64, 4, 4)
        model_custom = NiPiZD.instantiate(
            bgc_type; palatability_matrix=correct, assimilation_efficiency_matrix=correct
        )
        @test required_biogeochemical_tracers(model_custom) == (:N, :D, :Z1, :Z2, :P1, :P2)
    end

    @testset "Diameter overrides" begin
        @test_throws ArgumentError NiPiZD.instantiate(bgc_type; phyto_diameters=[1.0])
        @test_throws ArgumentError NiPiZD.instantiate(
            bgc_type; phyto_diameters=[1.0, 2.0, 3.0]
        )

        bgc_1p = NiPiZD.construct(; n_phyto=1)
        model_1p = NiPiZD.instantiate(bgc_1p; phyto_diameters=[2.0])
        @test required_biogeochemical_tracers(model_1p) == (:N, :D, :Z1, :Z2, :P1)
    end

    @testset "Geider-style growth" begin
        bgc_geider = NiPiZD.construct(;
            phyto_pft_parameters=NiPiZD.default_phyto_geider_pft_parameters(Float64),
            nutrient_dynamics=nutrient_geider_light,
            phyto_dynamics=phytoplankton_geider_light,
        )

        model_geider = NiPiZD.instantiate(bgc_geider)
        ordered_geider = [
            tracer_vals(s) for s in required_biogeochemical_tracers(model_geider)
        ]
        @test isfinite(model_geider(Val(:N), 0, 0, 0, 0, ordered_geider..., PAR))
    end

    @testset "Parameter variation" begin
        p0 = NiPiZD.default_phyto_pft_parameters(Float64)

        function with_growth_rate_a(p, new_a)
            return PhytoPFTParameters{Float64}(
                new_a,
                p.maximum_growth_rate_b,
                p.nutrient_half_saturation_a,
                p.nutrient_half_saturation_b,
                p.linear_mortality,
                p.alpha,
                p.photosynthetic_slope,
                p.chlorophyll_to_carbon_ratio,
                p.can_eat,
                p.can_be_eaten,
                p.optimum_predator_prey_ratio,
                p.protection,
                p.specificity,
                p.assimilation_efficiency,
            )
        end

        bgc_a5 = NiPiZD.construct(; phyto_pft_parameters=with_growth_rate_a(p0, 5 / day))
        bgc_a10 = NiPiZD.construct(; phyto_pft_parameters=with_growth_rate_a(p0, 10 / day))

        m5 = bgc_a5()
        m10 = bgc_a10()

        ordered5 = [tracer_vals(s) for s in required_biogeochemical_tracers(m5)]
        ordered10 = [tracer_vals(s) for s in required_biogeochemical_tracers(m10)]

        @test m5(Val(:P1), 0, 0, 0, 0, ordered5..., PAR) !=
            m10(Val(:P1), 0, 0, 0, 0, ordered10..., PAR)
    end

    @testset "Tracer sinking" begin
        @test_throws ArgumentError NiPiZD.instantiate(
            bgc_type; sinking_tracers=(P1=0.2551 / day, P2=0.2551 / day, D=2.7489 / day)
        )

        bgc_sink = NiPiZD.construct(;
            sinking_tracers=(P1=0.2551 / day, P2=0.2551 / day, D=2.7489 / day)
        )
        model_sink = bgc_sink()

        @test biogeochemical_drift_velocity(model_sink, Val(:P1)).w.data[1, 1, 1] ==
            -0.2551 / day
        @test biogeochemical_drift_velocity(model_sink, Val(:P2)).w.data[1, 1, 1] ==
            -0.2551 / day
        @test biogeochemical_drift_velocity(model_sink, Val(:D)).w.data[1, 1, 1] ==
            -2.7489 / day
        @test biogeochemical_drift_velocity(model_sink, Val(:Z1)).w == ZeroField()
    end
end
