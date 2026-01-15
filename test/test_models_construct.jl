using Agate
using Test

using Agate.Constructor: construct
using Agate.Models: NiPiZDFactory, DarwinFactory

using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    biogeochemical_drift_velocity


@testset "Agate.Constructor.construct" begin
    @testset "NiPiZD defaults" begin
        bgc = construct(NiPiZDFactory(); FT=Float32)

        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1, :P2)

        P1 = 0.01f0
        P2 = 0.01f0
        Z1 = 0.05f0
        Z2 = 0.05f0
        N  = 7.0f0
        D  = 1.0f0
        PAR = 100.0f0

        tracer_vals(sym) = sym === :P1 ? P1 :
                           sym === :P2 ? P2 :
                           sym === :Z1 ? Z1 :
                           sym === :Z2 ? Z2 :
                           sym === :N  ? N  :
                           D

        ordered = [tracer_vals(s) for s in required_biogeochemical_tracers(bgc)]

        @test isfinite(bgc(Val(:N),  0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:D),  0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:P1), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:Z1), 0, 0, 0, 0, ordered..., PAR))
    end

    @testset "NiPiZD interaction overrides" begin
        bgc = construct(NiPiZDFactory(); FT=Float32)
        n = size(bgc.parameters.palatability_matrix, 1)

        wrong = zeros(Float32, 2, 2)
        @test_throws ArgumentError construct(
            NiPiZDFactory();
            FT=Float32,
            interactions=(; palatability_matrix=wrong, assimilation_matrix=wrong),
        )

        correct = zeros(Float32, n, n)
        bgc2 = construct(
            NiPiZDFactory();
            FT=Float32,
            interactions=(; palatability_matrix=correct, assimilation_matrix=correct),
        )
        @test required_biogeochemical_tracers(bgc2) == (:N, :D, :Z1, :Z2, :P1, :P2)
    end

    @testset "NiPiZD community override" begin
        factory = NiPiZDFactory()
        base = Agate.Models.default_community(factory)

        # Override phytoplankton to a single explicit diameter.
        P = (; base.P..., diameters=[3.0], n=1)
        community = (Z = base.Z, P = P)

        bgc = construct(factory; community=community)
        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1)
    end

    @testset "NiPiZD sinking" begin
        bgc = construct(
            NiPiZDFactory();
            sinking_tracers=(P1=0.2551 / day, P2=0.2551 / day, D=2.7489 / day),
        )

        @test biogeochemical_drift_velocity(bgc, Val(:P1)).w.data[1, 1, 1] == -0.2551 / day
        @test biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1] == -2.7489 / day
        @test biogeochemical_drift_velocity(bgc, Val(:Z1)).w == ZeroField()
    end

    @testset "DARWIN defaults" begin
        bgc = construct(DarwinFactory(); FT=Float32)

        @test required_biogeochemical_tracers(bgc)[1:9] ==
              (:DIC, :DIN, :PO4, :DOC, :POC, :DON, :PON, :DOP, :POP)
    end

    @testset "GPU smoke test" begin
        # NOTE: Loading CUDA can crash Julia in misconfigured environments (e.g. mixed system/toolkit libs).
        # To keep the default test suite robust, this test only runs when explicitly enabled.
        if lowercase(get(ENV, "AGATE_TEST_CUDA", "0")) in ("1", "true", "yes")
            @eval using CUDA
            @eval using Oceananigans.Architectures: GPU, array_type

            if CUDA.functional()
                bgc_cpu = construct(NiPiZDFactory(); FT=Float32)
                bgc_gpu = construct(NiPiZDFactory(); FT=Float32, arch=GPU())

                @test required_biogeochemical_tracers(bgc_gpu) == required_biogeochemical_tracers(bgc_cpu)
                @test bgc_gpu.parameters.palatability_matrix isa array_type(GPU())
                @test bgc_gpu.parameters.maximum_predation_rate isa array_type(GPU())
            else
                @test true
            end
        else
            @test true
        end
    end

    @testset "Input validation" begin
        factory = NiPiZDFactory()

        base_args = Agate.Models.default_community(factory)
        base_dyn  = Agate.Models.default_plankton_dynamics(factory)

        # Missing group in community should produce a single ArgumentError.
        bad_args = (P = base_args.P,)  # missing Z
        @test_throws ArgumentError construct(factory; plankton_dynamics=base_dyn, community=bad_args)

        # Diameter range without `n` should error.
        bad_args2 = (Z = base_args.Z, P = (; base_args.P..., diameters=(2.0, 10.0, :log_splitting), n=nothing))
        @test_throws ArgumentError construct(factory; community=bad_args2)

        # Interaction matrix wrong size should error.
        bgc = construct(factory)
        n_total = size(bgc.parameters.palatability_matrix, 1)
        wrong = zeros(Float64, n_total + 1, n_total + 1)
        interactions = (; palatability_matrix=wrong, assimilation_matrix=wrong)
        @test_throws ArgumentError construct(factory; interactions=interactions)
        @test_throws ArgumentError construct(factory; interactions=(; palatability_matrix=wrong, assimilation_matrix=wrong))
    end
end
