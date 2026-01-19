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

        wrong = zeros(Float32, 3, 3)
        @test_throws ArgumentError NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=wrong,
            assimilation_matrix=wrong,
        )

        # Group-block matrices are expanded to full (n_total, n_total) matrices.
        # For role-aware interaction matrices, (n_groups, n_groups) matrices are
        # ambiguous, so we wrap them explicitly.
        block = Agate.Utils.GroupBlockMatrix(Float32[0 1; 2 3])
        bgc_block = NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=block,
            assimilation_matrix=block,
        )
        @test size(bgc_block.parameters.palatability_matrix) == (n, n)
        @test bgc_block.parameters.palatability_matrix[1, 1] == 0f0
        @test bgc_block.parameters.palatability_matrix[1, 3] == 1f0
        @test bgc_block.parameters.palatability_matrix[3, 1] == 2f0
        @test bgc_block.parameters.palatability_matrix[end, end] == 3f0

        # Rectangular consumer-by-prey matrices (n_consumer, n_prey) are embedded
        # into the full (n_total, n_total) storage with zeros elsewhere.
        rect = reshape(Float32.(1:4), 2, 2)
        bgc_rect = NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=rect,
            assimilation_matrix=rect,
        )
        M = bgc_rect.parameters.palatability_matrix
        @test size(M) == (n, n)
        @test M[1, 3] == 1f0
        @test M[2, 3] == 2f0
        @test M[1, 4] == 3f0
        @test M[2, 4] == 4f0
        @test all(M[3:4, :] .== 0f0)
        @test all(M[:, 1:2] .== 0f0)

        # Axis-local group-block matrices (n_consumer_groups, n_prey_groups) are
        # expanded within the consumer/prey axes and then embedded.
        axis_block = reshape(Float32[7], 1, 1)
        bgc_axis_block = NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=axis_block,
            assimilation_matrix=axis_block,
        )
        M2 = bgc_axis_block.parameters.palatability_matrix
        @test all(M2[1:2, 3:4] .== 7f0)
        @test all(M2[3:4, :] .== 0f0)
        @test all(M2[:, 1:2] .== 0f0)

        # Providers can return axis-sized rectangular matrices.
        rect_provider(ctx) = fill(Float32(9), length(ctx.consumer_indices), length(ctx.prey_indices))
        bgc_rect_provider = NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=rect_provider,
            assimilation_matrix=rect_provider,
        )
        M3 = bgc_rect_provider.parameters.palatability_matrix
        @test all(M3[1:2, 3:4] .== 9f0)
        @test all(M3[3:4, :] .== 0f0)
        @test all(M3[:, 1:2] .== 0f0)

        correct = zeros(Float32, n, n)
        bgc2 = NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=correct,
            assimilation_matrix=correct,
        )
        @test required_biogeochemical_tracers(bgc2) == (:N, :D, :Z1, :Z2, :P1, :P2)

        # Providers are allowed for the public keywords. They are evaluated once during construction.
        pal_provider(ctx) = zeros(Float32, ctx.n_total, ctx.n_total)
        assim_provider(ctx) = zeros(Float32, ctx.n_total, ctx.n_total)

        # Old-style providers are rejected.
        old_style_provider(diameters, group_symbols) = zeros(Float32, length(diameters), length(diameters))
        @test_throws ArgumentError NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=pal_provider,
            assimilation_matrix=old_style_provider,
        )

        bgc3 = NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=pal_provider,
            assimilation_matrix=assim_provider,
        )
        @test size(bgc3.parameters.palatability_matrix) == (n, n)
        @test size(bgc3.parameters.assimilation_matrix) == (n, n)


        # Providers may also return consumer-by-prey rectangular matrices when the
        # parameter directory declares axes for the matrix.
        rect_provider(ctx) = fill(Float32(9), length(ctx.consumer_indices), length(ctx.prey_indices))
        bgc_rect_provider = NiPiZD.construct(
            ;
            grid=dummy_grid(Float32),
            palatability_matrix=rect_provider,
            assimilation_matrix=rect_provider,
        )
        M3 = bgc_rect_provider.parameters.palatability_matrix
        @test all(M3[1:2, 3:4] .== 9f0)
        @test all(M3[3:4, :] .== 0f0)
        @test all(M3[:, 1:2] .== 0f0)
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
