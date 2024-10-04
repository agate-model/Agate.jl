using Agate.Models.Dynamic: expression_check
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry,
                                    required_biogeochemical_tracers,
                                    required_biogeochemical_auxiliary_fields

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
            parameters = (α=2/3, β=4/3, δ=1, γ=1)
            name = create_bgc_struct(:name, parameters)

            @test typeof(name) == DataType
            @test fieldnames(name) == (:α, :β, :δ, :γ)
        end
    end

    @testset "add_bgc_methods" begin

        @testset "all methods exist and behave as expected" begin

            parameters = (α=2/3, β=4/3, δ=1, γ=1)
            tracers = Dict(
                "R" => :(α*R - β*R*F),
                "F" => :(-γ*F + δ*R*F)
            )
            auxiliary_fields = [:PAR,]

            LV = create_bgc_struct(:LV, parameters)
            add_bgc_methods(LV, tracers, auxiliary_fields=auxiliary_fields)

            # instantiate the same model with different parameters
            model1 = LV()
            model2 = LV(1, 1, 2, 2)

            @test all(required_biogeochemical_tracers(model1) .=== [:R, :F])
            @test all(required_biogeochemical_auxiliary_fields(model1) .=== [:PAR])

            # the inputs to the tracer methods are: x,y,z,t,R,F,PAR
            @test model1(Val(:R), 0, 0, 0, 0, 10, 2, 0) == 2/3 * 10 - 4/3 * 10 * 2
            @test model1(Val(:F), 0, 0, 0, 0, 10, 2, 0) == -1 * 2 + 1 * 10 * 2

            @test model2(Val(:R), 0, 0, 0, 0, 10, 2, 0) == 1 * 10 - 1 * 10 * 2
            @test model2(Val(:F), 0, 0, 0, 0, 10, 2, 0) == -2 * 2 + 2 * 10 * 2

        end

        @testset "use helper functions" begin

            # NPZD model

            helper_functions = "./fixtures/NPZD.jl"

            parameters = (
                μ₀ = 0.6989,
                kₙ = 2.3868,
                lᵖⁿ = 0.066,
                lᶻⁿ = 0.0102,
                lᵖᵈ = 0.0101,
                gₘₐₓ = 2.1522,
                kₚ = 0.5573,
                β = 0.9116,
                lᶻᵈ = 0.3395,
                rᵈⁿ = 0.1213,
                α = 0.1953,
            )
            aux_field_vars = [:PAR]

            tracers = Dict(
                "N" => :(phytoplankton_metabolic_loss(P, lᵖⁿ)
                + zooplankton_metabolic_loss(Z, lᶻⁿ)
                + remineralization(D, rᵈⁿ)
                - phytoplankton_growth(N, P, PAR, μ₀, kₙ, α)),

                "D" => :(phytoplankton_mortality_loss(P, lᵖᵈ)
                + zooplankton_assimilation_loss(P, Z, β, gₘₐₓ, kₚ)
                + zooplankton_mortality_loss(Z, lᶻᵈ)
                - remineralization(D, rᵈⁿ)),

                "P" => :(phytoplankton_growth(N, P, PAR, μ₀, kₙ, α)
                - phytoplankton_grazing_loss(P, Z, gₘₐₓ, kₚ)
                - phytoplankton_metabolic_loss(P, lᵖⁿ)
                - phytoplankton_mortality_loss(P, lᵖᵈ)),

                "Z" => :(zooplankton_grazing_gain(P, Z, β, gₘₐₓ, kₚ)
                - zooplankton_metabolic_loss(Z, lᶻⁿ)
                - zooplankton_mortality_loss(Z, lᶻᵈ))
            )

            NPZD = create_bgc_struct(:NPZD, parameters)
            add_bgc_methods(NPZD, tracers, auxiliary_fields=aux_field_vars, helper_functions=helper_functions)
            model = NPZD()

            Z = 0.05
            P = 0.01
            N = 7.0
            D = 1
            PAR = 1

            @test isapprox(model(Val(:N), 0, 0, 0, 0, Z,P,N,D,PAR), 0.12106732850635774)
            @test isapprox(model(Val(:D), 0, 0, 0, 0, Z,P,N,D,PAR), -0.12034718812821411)
            @test isapprox(model(Val(:P), 0, 0, 0, 0, Z,P,N,D,PAR), 0.000607034934978301)
            @test isapprox(model(Val(:Z), 0, 0, 0, 0, Z,P,N,D,PAR), -0.0013271753131219258)

        end

    end

end
