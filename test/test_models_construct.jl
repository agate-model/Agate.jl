using Agate
using Test

using OceanBioME: BoxModelGrid

using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    biogeochemical_drift_velocity


@testset "Public model constructors" begin
    @testset "NiPiZD defaults" begin
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))

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
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))
        n = size(bgc.parameters.palatability_matrix, 1)

        wrong = zeros(Float32, 2, 2)
        @test_throws ArgumentError NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=wrong,
            assimilation_matrix=wrong,
        )

        correct = zeros(Float32, n, n)
        bgc2 = NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=correct,
            assimilation_matrix=correct,
        )
        @test required_biogeochemical_tracers(bgc2) == (:N, :D, :Z1, :Z2, :P1, :P2)

        # Bare matrix functions are intentionally unsupported: wrap in MatrixFn(f; deps=...).
        bare = (ctx, depvals) -> zeros(Float32, n, n)
        err = try
            NiPiZD.construct(
                ;
                grid=dummy_grid(Float32),
                interactions=(; palatability_matrix=bare, assimilation_matrix=correct),
            )
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("MatrixFn", sprint(showerror, err))
    end

    @testset "NiPiZD community structure overrides" begin
        bgc = NiPiZD.construct(; n_phyto=1, phyto_diameters=[3.0], grid=dummy_grid(Float32))
        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1)
    end

    @testset "NiPiZD sinking" begin
        bgc = NiPiZD.construct(
            ;
            sinking_tracers=(P1=0.2551 / day, P2=0.2551 / day, D=2.7489 / day),
        )

        @test biogeochemical_drift_velocity(bgc, Val(:P1)).w.data[1, 1, 1] == -0.2551 / day
        @test biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1] == -2.7489 / day
        @test biogeochemical_drift_velocity(bgc, Val(:Z1)).w == ZeroField()
    end

    @testset "DARWIN defaults" begin
        bgc = DARWIN.construct(; grid=dummy_grid(Float32))

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
                bgc_cpu = NiPiZD.construct(; grid=dummy_grid(Float32))
                bgc_gpu = NiPiZD.construct(; grid=dummy_grid(Float32; arch=GPU()))

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
        @test_throws ArgumentError NiPiZD.construct(; n_phyto=0)
        @test_throws ArgumentError NiPiZD.construct(; n_zoo=0)

        # Grid determines precision.
        bgc_f32 = NiPiZD.construct(; grid=dummy_grid(Float32))
        @test bgc_f32.parameters.detritus_remineralization isa Float32

        # Wrong interaction matrix sizes should error.
        bgc = NiPiZD.construct(; grid=dummy_grid(Float64))
        n_total = size(bgc.parameters.palatability_matrix, 1)
        wrong = zeros(Float64, n_total + 1, n_total + 1)
        @test_throws ArgumentError NiPiZD.construct(; grid=dummy_grid(Float64), palatability_matrix=wrong, assimilation_matrix=wrong)
    end
end
