using Test

include("production_oracle_helpers.jl")

@testset "Production NiPiZD trajectory" begin
    reference_path = joinpath(@__DIR__, "reference", "nipizd_production_v0.10.1.csv")
    if !isfile(reference_path)
        @test_skip isfile(reference_path)
    else
        box_model, _ = build_nipizd_production_box(Float64)
        set_nipizd_production_state!(box_model, Float64)
        actual = production_trajectory(box_model, NIPIZD_PRODUCTION_TRACERS)
        reference = read_trajectory_reference(reference_path, NIPIZD_PRODUCTION_TRACERS)

        @test length(actual) == length(reference.rows)
        rtol = get(reference.metadata, :trajectory_rtol, 0.0)
        atol = get(reference.metadata, :trajectory_atol, 0.0)
        for i in eachindex(actual, reference.rows)
            @test actual[i].step == reference.rows[i].step
            @test actual[i].time == reference.rows[i].time
            for tracer in NIPIZD_PRODUCTION_TRACERS
                @test isapprox(
                    Float64(getproperty(actual[i], tracer)),
                    getproperty(reference.rows[i], tracer);
                    rtol,
                    atol,
                )
            end
        end
    end
end
