using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans.Units

const year = years = 365day

reference_path = joinpath(@__DIR__, "reference", "nipizd_box_model.csv")
reference_P = [
    parse(Float64, split(row, ',')[2]) for row in readlines(reference_path)
    if !isempty(row) && !startswith(row, '#') && row != "day,P"
]

@testset "box_model" begin
    @testset "NPZD box model" begin
        bgc_model = Biogeochemistry(
            AgateNPZD(parameters); light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
        )
        box_model = BoxModel(; biogeochemistry=bgc_model)
        set!(box_model; N=7.0, P=0.01, Z=0.05, D=0.0)

        Δt = 1day
        for i in 1:1000
            time_step!(box_model, Δt)
            if i % 10 == 0
                P = box_model.fields.P.data[1, 1, 1]
                @test P ≈ reference_P[i ÷ 10] rtol=1e-12 atol=0
            end
        end

        total_N = sum(name -> box_model.fields[name].data[1, 1, 1], (:N, :P, :Z, :D))
        @test total_N ≈ 7.06 rtol=1e-12 atol=1e-12
    end
end
