using Agate
using Test

using OceanBioME: BoxModelGrid

using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers, biogeochemical_drift_velocity

@testset "Public model constructors" begin
    @testset "NiPiZD defaults" begin
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))

        # Guardrail for GPU compilation: tracer callables must be concretely typed.
        @test !any(t -> t === Any, fieldtypes(typeof(bgc.tracer_functions)))

        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1, :P2)

        P1 = 0.01f0
        P2 = 0.01f0
        Z1 = 0.05f0
        Z2 = 0.05f0
        N = 7.0f0
        D = 1.0f0
        PAR = 100.0f0

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

        ordered = [tracer_vals(s) for s in required_biogeochemical_tracers(bgc)]

        @test isfinite(bgc(Val(:N), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:D), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:P1), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(bgc(Val(:Z1), 0, 0, 0, 0, ordered..., PAR))
    end

    @testset "NiPiZD interaction overrides" begin
        bgc = NiPiZD.construct(; grid=dummy_grid(Float32))
        ints0 = bgc.parameters.interactions
        n_total = length(ints0.global_to_prey)
        n_cons = length(ints0.consumer_global)
        n_prey = length(ints0.prey_global)

        wrong = zeros(Float32, 3, 3)
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32), palatability_matrix=wrong, assimilation_matrix=wrong
        )

        # Group-block matrices are accepted when wrapped explicitly. They are expanded
        # over all groups and then sliced to the consumer-by-prey axes.
        block = Agate.Configuration.GroupBlockMatrix(Float32[0 1; 2 3])
        bgc_block = NiPiZD.construct(;
            grid=dummy_grid(Float32), palatability_matrix=block, assimilation_matrix=block
        )
        M0 = bgc_block.parameters.palatability_matrix
        @test size(M0) == (n_cons, n_prey)
        @test all(M0 .== 1.0f0)

        # Rectangular consumer-by-prey matrices are stored as-is.
        rect = reshape(Float32.(1:(n_cons * n_prey)), n_cons, n_prey)
        bgc_rect = NiPiZD.construct(;
            grid=dummy_grid(Float32), palatability_matrix=rect, assimilation_matrix=rect
        )
        M = bgc_rect.parameters.palatability_matrix
        @test size(M) == (n_cons, n_prey)
        @test all(M .== rect)

        # Axis-local group-block matrices (consumer_groups, prey_groups) are expanded
        # within the consumer/prey axes.
        axis_block = reshape(Float32[7], 1, 1)
        bgc_axis_block = NiPiZD.construct(;
            grid=dummy_grid(Float32),
            palatability_matrix=axis_block,
            assimilation_matrix=axis_block,
        )
        M2 = bgc_axis_block.parameters.palatability_matrix
        @test size(M2) == (n_cons, n_prey)
        @test all(M2 .== 7.0f0)

        # Providers can return axis-sized rectangular matrices.
        rect_provider(ctx) =
            fill(Float32(9), length(ctx.consumer_indices), length(ctx.prey_indices))
        bgc_rect_provider = NiPiZD.construct(;
            grid=dummy_grid(Float32),
            palatability_matrix=rect_provider,
            assimilation_matrix=rect_provider,
        )
        M3 = bgc_rect_provider.parameters.palatability_matrix
        @test size(M3) == (n_cons, n_prey)
        @test all(M3 .== 9.0f0)

        # Full-square matrices are accepted but are stored sliced to the consumer/prey axes.
        correct = zeros(Float32, n_total, n_total)
        bgc2 = NiPiZD.construct(;
            grid=dummy_grid(Float32),
            palatability_matrix=correct,
            assimilation_matrix=correct,
        )
        @test required_biogeochemical_tracers(bgc2) == (:N, :D, :Z1, :Z2, :P1, :P2)
        @test size(bgc2.parameters.palatability_matrix) == (n_cons, n_prey)

        # Providers are allowed for the public keywords. They are evaluated once during construction.
        pal_provider(ctx) = zeros(Float32, ctx.n_total, ctx.n_total)
        assim_provider(ctx) = zeros(Float32, ctx.n_total, ctx.n_total)

        # Old-style providers are rejected.
        old_style_provider(diameters, group_symbols) =
            zeros(Float32, length(diameters), length(diameters))
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float32),
            palatability_matrix=pal_provider,
            assimilation_matrix=old_style_provider,
        )

        bgc3 = NiPiZD.construct(;
            grid=dummy_grid(Float32),
            palatability_matrix=pal_provider,
            assimilation_matrix=assim_provider,
        )
        @test size(bgc3.parameters.palatability_matrix) == (n_cons, n_prey)
        @test size(bgc3.parameters.assimilation_matrix) == (n_cons, n_prey)
    end

    @testset "Derived interaction matrices" begin
        # If a model exposes interaction traits, overriding one of those traits
        # should regenerate the derived matrices (unless the matrix itself is
        # explicitly overridden).

        bgc0 = NiPiZD.construct(; grid=dummy_grid(Float32))
        pal0 = bgc0.parameters.interactions.palatability
        n_total = length(bgc0.parameters.interactions.global_to_prey)

        specificity = zeros(Float32, n_total)
        specificity[bgc0.parameters.interactions.consumer_global] .= 3.0f0

        bgc1 = NiPiZD.construct(;
            grid=dummy_grid(Float32), parameters=(; specificity=specificity)
        )
        pal1 = bgc1.parameters.interactions.palatability
        @test any(pal1 .!= pal0)

        rect = fill(Float32(11), size(pal0))
        bgc2 = NiPiZD.construct(;
            grid=dummy_grid(Float32),
            parameters=(; specificity=specificity),
            palatability_matrix=rect,
        )
        @test all(bgc2.parameters.interactions.palatability .== rect)

        dar0 = DARWIN.construct(; grid=dummy_grid(Float32))
        dar_pal0 = dar0.parameters.interactions.palatability
        dar_n_total = length(dar0.parameters.interactions.global_to_prey)
        dar_spec = zeros(Float32, dar_n_total)
        dar_spec[dar0.parameters.interactions.consumer_global] .= 2.0f0
        dar1 = DARWIN.construct(;
            grid=dummy_grid(Float32), parameters=(; specificity=dar_spec)
        )
        dar_pal1 = dar1.parameters.interactions.palatability
        @test any(dar_pal1 .!= dar_pal0)
    end

    @testset "NiPiZD community structure overrides" begin
        bgc = NiPiZD.construct(; n_phyto=1, phyto_diameters=[3.0], grid=dummy_grid(Float32))
        @test required_biogeochemical_tracers(bgc) == (:N, :D, :Z1, :Z2, :P1)
    end

    @testset "NiPiZD sinking" begin
        bgc = NiPiZD.construct(;
            sinking_tracers=(P1=0.2551 / day, P2=0.2551 / day, D=2.7489 / day)
        )

        @test biogeochemical_drift_velocity(bgc, Val(:P1)).w.data[1, 1, 1] == -0.2551 / day
        @test biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1] == -2.7489 / day
        @test biogeochemical_drift_velocity(bgc, Val(:Z1)).w == ZeroField()
    end

    @testset "DARWIN defaults" begin
        bgc = DARWIN.construct(; grid=dummy_grid(Float32))

        @test !any(t -> t === Any, fieldtypes(typeof(bgc.tracer_functions)))

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

                @test required_biogeochemical_tracers(bgc_gpu) ==
                    required_biogeochemical_tracers(bgc_cpu)
                @test bgc_gpu.parameters.interactions.palatability isa array_type(GPU())
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
        n_cons = size(bgc.parameters.palatability_matrix, 1)
        n_prey = size(bgc.parameters.palatability_matrix, 2)
        wrong = zeros(Float64, n_cons + 1, n_prey + 1)
        @test_throws ArgumentError NiPiZD.construct(;
            grid=dummy_grid(Float64), palatability_matrix=wrong, assimilation_matrix=wrong
        )
    end

    @testset "InteractionBlocks helpers" begin
        using Agate.Configuration:
            roles_from_groups, interaction_blocks, set_block!, scale_block!, forbid_link!

        roles = roles_from_groups(; consumers=:Z, prey=(:P, :Z))

        pal = interaction_blocks(roles; init=0)
        set_block!(pal; consumer_group=:Z, prey_group=:P, value=1.0f0)
        set_block!(pal; consumer_group=:Z, prey_group=:Z, value=0.25f0)
        forbid_link!(pal; consumer_group=:Z, prey_group=:Z)

        factory = Agate.NiPiZD.NiPiZDFactory()
        base = Agate.Factories.default_community(factory)
        community = Agate.Configuration.build_plankton_community(
            base;
            n=(Z=2, P=2),
            diameters=(Z=(20, 100, :linear_splitting), P=(2, 10, :log_splitting)),
        )

        bgc = Agate.Construction.construct_factory(
            factory;
            grid=dummy_grid(Float32),
            community=community,
            interaction_roles=roles,
            interaction_overrides=(palatability_matrix=pal,),
            auxiliary_fields=(:PAR,),
        )
        M = bgc.parameters.palatability_matrix
        @test size(M, 1) == 2
        @test size(M, 2) == 4
        @test all(M[:, 1:2] .== 0.0f0) # Z as prey
        @test all(M[:, 3:4] .== 1.0f0) # P as prey

        scale_block!(pal; consumer_group=:Z, prey_group=:P, factor=0.5f0)
        bgc2 = Agate.Construction.construct_factory(
            factory;
            grid=dummy_grid(Float32),
            community=community,
            interaction_roles=roles,
            interaction_overrides=(palatability_matrix=pal,),
            auxiliary_fields=(:PAR,),
        )
        M2 = bgc2.parameters.palatability_matrix
        @test all(M2[:, 3:4] .== 0.5f0)
    end
end
