using Agate
using Agate.Models.Dynamic: expression_check, create_bgc_struct, add_bgc_methods!
using OceanBioME
using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    AbstractContinuousFormBiogeochemistry,
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

@testset "Models.Dynamic" begin
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

        # no errors - method defined in Dynamic module
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
            parameters = (α=2 / 3, β=4 / 3, δ=1, γ=1)
            name = create_bgc_struct(:name, parameters)

            @test typeof(name) == DataType
            @test fieldnames(name) == (:α, :β, :δ, :γ)
        end
    end

    @testset "add_bgc_methods" begin
        @testset "core methods exist and behave as expected" begin
            parameters = (α=2 / 3, β=4 / 3, δ=1, γ=1)
            tracers = Dict("R" => :(α * R - β * R * F), "F" => :(-γ * F + δ * R * F))
            auxiliary_fields = [:PAR]

            LV = create_bgc_struct(:LV, parameters)
            add_bgc_methods!(LV, tracers; auxiliary_fields=auxiliary_fields)

            # instantiate the same model with different parameters
            model1 = LV()
            model2 = LV(1, 1, 2, 2)

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
            include(joinpath("..", "examples", "NPZD", "tracers.jl"))
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
            include(joinpath("..", "examples", "NPZD", "tracers.jl"))

            NPZD_sink = define_tracer_functions(
                parameters,
                tracers;
                helper_functions=joinpath("..", "examples", "NPZD", "functions.jl"),
                sinking_tracers=(P=0.2551 / day, D=2.7489 / day),
                grid=BoxModelGrid(),
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
