using Agate
using Test

using Agate.Models: construct, NiPiZDFactory, DarwinFactory
using Agate.Utils: pft_get

using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    biogeochemical_drift_velocity

using Adapt

@testset "Agate.Models.construct" begin
    @testset "NiPiZD defaults" begin
        bgc_type = construct(NiPiZDFactory(); FT=Float32)
        bgc = bgc_type()

        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1, :P2)

        P1 = 0.01f0
        P2 = 0.01f0
        Z1 = 0.05f0
        Z2 = 0.05f0
        N = 7.0f0
        D = 1.0f0
        PAR = 100.0f0

        tracer_vals(sym) = sym === :P1 ? P1 :
                           sym === :P2 ? P2 :
                           sym === :Z1 ? Z1 :
                           sym === :Z2 ? Z2 :
                           sym === :N ? N :
                           D

        ordered = [tracer_vals(s) for s in required_biogeochemical_tracers(bgc)]

        @test isfinite(bgc(Val(:N), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:D), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:P1), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:Z1), 0, 0, 0, 0, ordered..., PAR))
    end

    @testset "NiPiZD interaction overrides" begin
        bgc_type = construct(NiPiZDFactory(); FT=Float32)
        bgc = bgc_type()
        n = length(bgc.parameters.data.diameters)

        wrong = zeros(Float32, 2, 2)
        @test_throws ArgumentError construct(
            NiPiZDFactory();
            FT=Float32,
            interactions=(; palatability_matrix=wrong, assimilation_efficiency_matrix=wrong),
        )

        correct = zeros(Float32, n, n)
        bgc_type2 = construct(
            NiPiZDFactory();
            FT=Float32,
            interactions=(; palatability_matrix=correct, assimilation_efficiency_matrix=correct),
        )
        @test required_biogeochemical_tracers(bgc_type2()) == (:N, :D, :Z1, :Z2, :P1, :P2)
    end

    @testset "NiPiZD plankton_args override" begin
        factory = NiPiZDFactory()
        base = Agate.Models.default_plankton_args(factory, Float64)
        # Override phytoplankton to a single explicit diameter.
        P = (; base.P..., diameters=[3.0], n=1)
        plankton_args = (Z = base.Z, P = P)

        bgc_type = construct(factory; FT=Float64, plankton_args=plankton_args)
        @test required_biogeochemical_tracers(bgc_type()) == (:N, :D, :Z1, :Z2, :P1)
    end

    @testset "NiPiZD sinking" begin
        bgc_type = construct(
            NiPiZDFactory();
            FT=Float64,
            sinking_tracers=(P1=0.2551 / day, P2=0.2551 / day, D=2.7489 / day),
        )

        bgc = bgc_type()
        @test biogeochemical_drift_velocity(bgc, Val(:P1)).w.data[1, 1, 1] == -0.2551 / day
        @test biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1] == -2.7489 / day
        @test biogeochemical_drift_velocity(bgc, Val(:Z1)).w == ZeroField()
    end

    @testset "DARWIN defaults" begin
        bgc_type = construct(DarwinFactory(); FT=Float32)
        bgc = bgc_type()

        @test required_biogeochemical_tracers(bgc)[1:9] ==
            (:DIC, :DIN, :PO4, :DOC, :POC, :DON, :PON, :DOP, :POP)
    end

    @testset "GPU adapt smoke test" begin
        # NOTE: Loading CUDA can crash Julia in misconfigured environments (e.g. mixed system/toolkit libs).
        # To keep the default test suite robust, this test only runs when explicitly enabled.
        if lowercase(get(ENV, "AGATE_TEST_CUDA", "0")) in ("1", "true", "yes")
            bgc_type = construct(NiPiZDFactory(); FT=Float32)
            bgc = bgc_type()

            @eval using CUDA
            @eval using Adapt

            if CUDA.functional()
                adapted = Adapt.adapt(CUDA.CuArray, bgc)
                @test required_biogeochemical_tracers(adapted) == required_biogeochemical_tracers(bgc)
            else
                @test true
            end
        else
            @test true
        end
    end

    @testset "Input validation" begin
        factory = NiPiZDFactory()

        # Missing group in plankton_args should produce a single ArgumentError.
        base_args = Agate.Models.default_plankton_args(factory, Float64)
        base_dyn  = Agate.Models.default_plankton_dynamics(factory)
        bad_args = (P = base_args.P,)  # missing Z
        @test_throws ArgumentError construct(factory; FT=Float64, plankton_dynamics=base_dyn, plankton_args=bad_args)

        # Diameter range without `n` should error.
        bad_P = (; base_args.P..., n=nothing)
        bad_args2 = (Z = base_args.Z, P = (; base_args.P..., diameters=(2.0, 10.0, :log_splitting)))
        # remove n from P spec entirely
        bad_args2 = (Z = base_args.Z, P = (; base_args.P..., diameters=(2.0, 10.0, :log_splitting), n=nothing))
        @test_throws ArgumentError construct(factory; FT=Float64, plankton_args=bad_args2)

        # Interaction matrix wrong size should error.
        bgc_type = construct(factory; FT=Float64)
        bgc = bgc_type()
        n_total = length(bgc.parameters.data.diameters)
        wrong = zeros(Float64, n_total + 1, n_total + 1)
        interactions = (; palatability_matrix=wrong, assimilation_efficiency_matrix=wrong)
        @test_throws ArgumentError construct(factory; FT=Float64, interactions=interactions)
        @test_throws ArgumentError construct(factory; FT=Float64, interactions=(; palatability_matrix=wrong, assimilation_efficiency_matrix=wrong))

    end
end
