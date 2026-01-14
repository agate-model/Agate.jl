using Agate
using Agate.Constructor:
    expression_check, create_bgc_struct, add_bgc_methods!, define_tracer_functions

using OceanBioME
using Oceananigans.Units
using Oceananigans.Fields: ZeroField

using OceanBioME: setup_velocity_fields
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

@testset "Utils" begin
    @testset "expression_check" begin
        f_expr = :(α * x - β * x * y)
        allowed = (:α, :β)
        @test_throws UndefVarError expression_check(allowed, f_expr)

        f_expr = :(α * x - β * x * y)
        allowed = (:α, :β, :x, :y)
        @test expression_check(allowed, f_expr) === nothing

        f_expr = :(create_bgc_struct(sn, p))
        allowed = (:sn, :p)
        @test expression_check(allowed, f_expr) === nothing
    end

    @testset "create_bgc_struct and add_bgc_methods" begin
        struct LVParameters{FT<:AbstractFloat}
            α::FT
            β::FT
            δ::FT
            γ::FT
        end

        parameters = LVParameters{Float64}(2 / 3, 4 / 3, 1.0, 1.0)
        tracers = (R=:(α * R - β * R * F), F=:(-γ * F + δ * R * F))
        auxiliary_fields = (:PAR,)

        LV = create_bgc_struct(:LV, parameters)
        add_bgc_methods!(LV, tracers, parameters; auxiliary_fields=auxiliary_fields)

        model1 = LV(parameters)
        model2 = LV(LVParameters{Float64}(1.0, 1.0, 2.0, 2.0))

        @test required_biogeochemical_tracers(model1) == (:R, :F)
        @test required_biogeochemical_auxiliary_fields(model1) == (:PAR,)

        # Method signature is: x, y, z, t, R, F, PAR
        @test model1(Val(:R), 0, 0, 0, 0, 10, 2, 0) == 2 / 3 * 10 - 4 / 3 * 10 * 2
        @test model1(Val(:F), 0, 0, 0, 0, 10, 2, 0) == -1 * 2 + 1 * 10 * 2

        @test model2(Val(:R), 0, 0, 0, 0, 10, 2, 0) == 1 * 10 - 1 * 10 * 2
        @test model2(Val(:F), 0, 0, 0, 0, 10, 2, 0) == -2 * 2 + 2 * 10 * 2
    end

    @testset "helper functions and tracer sinking" begin
        include(joinpath("NPZD", "tracers.jl"))

        model = AgateNPZD(parameters)

        Z = 0.05
        P = 0.01
        N = 7.0
        D = 1.0
        PAR = 1.0

        tracer_vals(sym) =
            if sym === :Z
                Z
            elseif sym === :P
                P
            elseif sym === :N
                N
            else
                D
            end
        ordered = [tracer_vals(s) for s in required_biogeochemical_tracers(model)]

        @test isfinite(model(Val(:N), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(model(Val(:D), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(model(Val(:P), 0, 0, 0, 0, ordered..., PAR))
        @test isfinite(model(Val(:Z), 0, 0, 0, 0, ordered..., PAR))

        sinking_velocities = setup_velocity_fields(
            (P=0.2551 / day, D=2.7489 / day), BoxModelGrid(), true
        )

        NPZD_sink = define_tracer_functions(
            parameters,
            tracers;
            helper_functions=joinpath(@__DIR__, "NPZD", "functions.jl"),
            sinking_velocities=sinking_velocities,
        )

        model_sink = NPZD_sink(parameters, sinking_velocities)

        @test biogeochemical_drift_velocity(model_sink, Val(:P)).w.data[1, 1, 1] ==
            -0.2551 / day
        @test biogeochemical_drift_velocity(model_sink, Val(:D)).w.data[1, 1, 1] ==
            -2.7489 / day
        @test biogeochemical_drift_velocity(model_sink, Val(:Z)).w == ZeroField()
    end
end
