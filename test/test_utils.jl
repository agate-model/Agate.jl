using Agate
using Test

using Agate.Functors: CompiledEquation, req
using Agate.Constructor: define_tracer_functions

using OceanBioME: BoxModelGrid, setup_velocity_fields
using Oceananigans.Units
using Oceananigans.Fields: ZeroField
using Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

@testset "Utils" begin
    @testset "define_tracer_functions" begin
        parameters = (α=2 / 3, β=4 / 3, δ=1.0, γ=1.0)

        fR = (bgc, x, y, z, t, R, F) -> begin
            p = bgc.parameters
            p.α * R - p.β * R * F
        end

        fF = (bgc, x, y, z, t, R, F) -> begin
            p = bgc.parameters
            -p.γ * F + p.δ * R * F
        end

        tracers = (
            R=CompiledEquation(fR, req(; scalars=(:α, :β))),
            F=CompiledEquation(fF, req(; scalars=(:γ, :δ))),
        )

        bgc_type = define_tracer_functions(parameters, tracers; auxiliary_fields=())
        bgc = bgc_type(parameters)

        @test required_biogeochemical_tracers(bgc) == (:R, :F)
        @test required_biogeochemical_auxiliary_fields(bgc) == ()

        R = 0.5
        F = 0.25

        @test isfinite(bgc(Val(:R), 0, 0, 0, 0, R, F))
        @test isfinite(bgc(Val(:F), 0, 0, 0, 0, R, F))
    end

    @testset "define_tracer_functions error context" begin
        parameters = (α=2.0,)
        fR = (bgc, x, y, z, t, R) -> begin
            p = bgc.parameters
            p.α * R
        end

        tracers = (R=CompiledEquation(fR, req(; scalars=(:α, :missing_param))),)

        err = try
            define_tracer_functions(parameters, tracers; auxiliary_fields=())
            nothing
        catch e
            e
        end

        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("Tracer :R", msg)
        @test occursin("missing_param", msg)
    end

    @testset "Optional sinking velocities" begin
        parameters = (k=1.0,)
        fC = (bgc, x, y, z, t, C) -> begin
            p = bgc.parameters
            -p.k * C
        end

        tracers = (C=CompiledEquation(fC, req(; scalars=(:k,))),)

        grid = BoxModelGrid()
        sinking_tracers = (C=1.0 / day,)
        sinking_velocities = setup_velocity_fields(sinking_tracers, grid, true)

        bgc_type = define_tracer_functions(
            parameters, tracers; auxiliary_fields=(), sinking_velocities=sinking_velocities
        )

        bgc = bgc_type(parameters, sinking_velocities)
        @test biogeochemical_drift_velocity(bgc, Val(:C)).w.data[1, 1, 1] == -1.0 / day
        @test biogeochemical_drift_velocity(bgc, Val(:missing_tracer)).w == ZeroField()
    end
end
