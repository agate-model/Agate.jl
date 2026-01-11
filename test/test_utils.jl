using Agate
using Test

using Agate.Constructor: ModelSpecification
using Agate.Utils: define_tracer_functions, expression_check

using OceanBioME: BoxModelGrid, setup_velocity_fields
using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

@testset "Utils" begin
    @testset "expression_check" begin
        f_expr = :(α * x - β * x * y)
        allowed = (:α, :β)
        @test_throws UndefVarError expression_check(allowed, f_expr)

        allowed = (:α, :β, :x, :y)
        @test expression_check(allowed, f_expr) === nothing
    end

    @testset "define_tracer_functions with ModelSpecification" begin
        parameters = ModelSpecification((α = 2 / 3, β = 4 / 3, δ = 1.0, γ = 1.0))
        tracers = (R = :(α * R - β * R * F), F = :(-γ * F + δ * R * F))

        bgc_type = define_tracer_functions(parameters, tracers; auxiliary_fields=())
        bgc = bgc_type()

        @test required_biogeochemical_tracers(bgc) == (:R, :F)
        @test required_biogeochemical_auxiliary_fields(bgc) == ()

        R = 0.5
        F = 0.25

        @test isfinite(bgc(Val(:R), 0, 0, 0, 0, R, F))
        @test isfinite(bgc(Val(:F), 0, 0, 0, 0, R, F))
    end

    @testset "Optional sinking velocities" begin
        parameters = ModelSpecification((k = 1.0,))
        tracers = (C = :(-k * C),)

        grid = BoxModelGrid()
        sinking_tracers = (C = 1.0 / day,)
        sinking_velocities = setup_velocity_fields(sinking_tracers, grid, true)

        bgc_type = define_tracer_functions(
            parameters,
            tracers;
            auxiliary_fields=(),
            sinking_velocities=sinking_velocities,
        )

        bgc = bgc_type()
        @test biogeochemical_drift_velocity(bgc, Val(:C)).w.data[1, 1, 1] == -1.0 / day
        @test biogeochemical_drift_velocity(bgc, Val(:missing_tracer)).w == ZeroField()
    end
end
