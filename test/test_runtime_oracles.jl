using Test
using Agate
using OceanBioME
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using Oceananigans.Fields: ZeroField
using Oceananigans.Units: day

@testset "Runtime oracles" begin
    @testset "NiPiZD sinking realization" begin
        for T in (Float64, Float32), open_bottom in (true, false)
            grid = OceanBioME.BoxModelGrid(T)
            authored = (D=T(2.5 / day), P_1=T(0.25 / day))
            bgc = Agate.Models.NiPiZD.construct(; grid, sinking_tracers=authored, open_bottom)
            converted = NamedTuple{keys(authored)}(Tuple(convert(T, v) for v in values(authored)))
            expected = OceanBioME.setup_velocity_fields(converted, grid, open_bottom)

            @test typeof(bgc.sinking_velocities) === typeof(expected)
            @test bgc.sinking_velocities.D.data == expected.D.data
            @test bgc.sinking_velocities.P_1.data == expected.P_1.data
            @test eltype(bgc.sinking_velocities.D.data) === T
            if hasproperty(expected.D, :boundary_conditions)
                @test repr(bgc.sinking_velocities.D.boundary_conditions) ==
                    repr(expected.D.boundary_conditions)
            end
            @test biogeochemical_drift_velocity(bgc, Val(:D)).w.data == expected.D.data
            @test biogeochemical_drift_velocity(bgc, Val(:P_1)).w.data == expected.P_1.data
            @test biogeochemical_drift_velocity(bgc, Val(:Z_1)).w == ZeroField()
        end

        bgc = Agate.Models.NiPiZD.construct(; sinking_tracers=(D=2.5 / day,))
        @test bgc.parameters.detritus_remineralization isa Float64
        @test biogeochemical_drift_velocity(bgc, Val(:D)).w.data[1, 1, 1] == -2.5 / day
    end

    @testset "Warmed local tendency" begin
        bgc = Agate.Models.NiPiZD.construct(; grid=dummy_grid(Float64))
        args = (0, 0, 0, 0, 7.0, 0.0, 0.05, 0.05, 0.01, 0.01, 100.0)
        tendency() = bgc(Val(:P_1), args...)

        tendency()
        tendency()
        value = @inferred tendency()
        allocations = @allocated tendency()

        @test isfinite(value)
        @test allocations == 0
    end
end
