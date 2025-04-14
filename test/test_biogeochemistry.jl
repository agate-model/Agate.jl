using Agate
using Agate.Utils: expression_check, create_bgc_struct, add_bgc_methods!

using OceanBioME
using Oceananigans.Units
using Oceananigans.Fields: ZeroField

using OceanBioME: setup_velocity_fields
using Oceananigans.Biogeochemistry:
    AbstractContinuousFormBiogeochemistry,
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

@testset "Utils" begin
    @testset "expression_check" begin

        # missing args
        f_expr = :(α * x - β * x * y)
        params = [:α, :β]
        @test_throws UndefVarError expression_check(params, f_expr)

        # method is not defined
        f_expr = :(f1(x, α) + f2(y, β))
        params = [:α, :β, :x, :y]
        @test_throws UndefVarError expression_check(params, f_expr)

        # no errors - using base methods
        f_expr = :(α * x - β * x * y)
        params = [:α, :β, :x, :y]
        @test expression_check(params, f_expr) === nothing

        # no errors - method defined Biogeochemistry module
        f_expr = :(create_bgc_struct(sn, p))
        params = [:sn, :p]
        @test expression_check(params, f_expr) === nothing

        # no errors - use of vectors
        f_expr = :(sum[a, b, c])
        params = [:a, :b, :c]
        @test expression_check(params, f_expr) === nothing
    end

    @testset "create_bgc_struct" begin
        @testset "invalid fieldnames throw an error" begin
            parameters = (x=1,)
            @test_throws DomainError create_bgc_struct(:name, parameters)

            parameters = (y=2,)
            @test_throws DomainError create_bgc_struct(:name, parameters)

            parameters = (z=3,)
            @test_throws DomainError create_bgc_struct(:name, parameters)

            parameters = (t=4,)
            @test_throws DomainError create_bgc_struct(:name, parameters)
        end

        @testset "data type created succesfully" begin
            parameters = (α=2 / 3, β=4 / 3, δ=1.0, γ=1.0)
            name = create_bgc_struct(:name, parameters)
            @test fieldnames(name) == (:α, :β, :δ, :γ)
        end
    end

    @testset "add_bgc_methods" begin
        @testset "core methods exist and behave as expected" begin
            parameters = (α=2 / 3, β=4 / 3, δ=1.0, γ=1.0)
            tracers = Dict("R" => :(α * R - β * R * F), "F" => :(-γ * F + δ * R * F))
            auxiliary_fields = [:PAR]

            LV = create_bgc_struct(:LV, parameters)
            add_bgc_methods!(LV, tracers; auxiliary_fields=auxiliary_fields)

            # instantiate the same model with different parameters
            model1 = LV()
            model2 = LV(1.0, 1.0, 2.0, 2.0)

            @test all(required_biogeochemical_tracers(model1) .=== [:R, :F])
            @test all(required_biogeochemical_auxiliary_fields(model1) .=== [:PAR])

            # the inputs to the tracer methods are: x,y,z,t,R,F,PAR
            @test model1(Val(:R), 0, 0, 0, 0, 10, 2, 0) == 2 / 3 * 10 - 4 / 3 * 10 * 2
            @test model1(Val(:F), 0, 0, 0, 0, 10, 2, 0) == -1 * 2 + 1 * 10 * 2

            @test model2(Val(:R), 0, 0, 0, 0, 10, 2, 0) == 1 * 10 - 1 * 10 * 2
            @test model2(Val(:F), 0, 0, 0, 0, 10, 2, 0) == -2 * 2 + 2 * 10 * 2
        end

        @testset "use helper functions" begin

            # NPZD model
            include(joinpath("NPZD", "tracers.jl"))
            model = NPZD()

            Z = 0.05
            P = 0.01
            N = 7.0
            D = 1
            PAR = 1

            @test isapprox(
                model(Val(:N), 0, 0, 0, 0, Z, P, N, D, PAR), 1.4012422280828442e-6
            )
            @test isapprox(
                model(Val(:D), 0, 0, 0, 0, Z, P, N, D, PAR), -1.3929072700024781e-6
            )
            @test isapprox(
                model(Val(:P), 0, 0, 0, 0, Z, P, N, D, PAR), 7.025867302989598e-9
            )
            @test isapprox(
                model(Val(:Z), 0, 0, 0, 0, Z, P, N, D, PAR), -1.5360825383355622e-8
            )
        end

        @testset "tracer sinking" begin
            include(joinpath("NPZD", "tracers.jl"))

            # if one uses BoxModelGrid and sets open_bottom to false then all velocities are
            # set to 0 (because the function smooths them to get to 0 when they reach the
            # bottom and in a BoxModel the tracers already are at "the bottom")
            sinking_velocities = setup_velocity_fields(
                (P=0.2551 / day, D=2.7489 / day), BoxModelGrid(), true
            )

            @show keys(sinking_velocities)
            @show sinking_velocities[:P]
            @show sinking_velocities[:P].data[1, 1, 1]

            NPZD_sink = define_tracer_functions(
                parameters,
                tracers;
                helper_functions=joinpath("NPZD", "functions.jl"),
                sinking_velocities=sinking_velocities,
            )

            model = NPZD_sink()

            @test OceanBioME.Models.Sediments.sinking_tracers(model) == (:P, :D)

            @test biogeochemical_drift_velocity(model, Val(:P)).w.data[1, 1, 1] ==
                -0.2551 / day
            @test biogeochemical_drift_velocity(model, Val(:D)).w.data[1, 1, 1] ==
                -2.7489 / day
            @test biogeochemical_drift_velocity(model, Val(:Z)).w == ZeroField()
        end
    end
end
