using Agate

const NiPiZD = Agate.Models.NiPiZD

using Agate.Library.Light
using Agate.Diagnostics: box_model_mass_balance

using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units: day, minutes

@testset "mass_balance" begin
    function build_box_model(bgc_instance, grid)
        bgc_model = Biogeochemistry(
            bgc_instance; light_attenuation=FunctionFieldPAR(; grid)
        )
        return BoxModel(; biogeochemistry=bgc_model, grid)
    end

    dt = 10minutes
    nsteps = Int(10day ÷ dt)

    @testset "NiPiZD box model" begin
        grid = BoxModelGrid()
        bgc_instance = NiPiZD.construct(; grid)
        box_model = build_box_model(bgc_instance, grid)
        set!(box_model; N=7, P_1=0.01, P_2=0.01, Z_1=0.05, Z_2=0.05, D=0.0)

        budgets = (total=[:N => 1, :P_1 => 1, :P_2 => 1, :Z_1 => 1, :Z_2 => 1, :D => 1],)

        result = box_model_mass_balance(box_model, budgets; dt, nsteps)
        @test isapprox(result.initial.total, result.final.total; rtol=1e-12, atol=0.0)
    end

    @testset "NiPiZD Float32 conservation" begin
        grid = BoxModelGrid(Float32)
        bgc_instance = NiPiZD.construct(; grid, scalar_type=Float32)
        box_model = build_box_model(bgc_instance, grid)
        set!(box_model; N=7f0, P_1=0.01f0, P_2=0.01f0, Z_1=0.05f0, Z_2=0.05f0, D=0f0)

        budgets = (total=[:N => 1, :P_1 => 1, :P_2 => 1, :Z_1 => 1, :Z_2 => 1, :D => 1],)

        result = box_model_mass_balance(box_model, budgets; dt, nsteps)
        @test result.initial.total isa Float32
        @test isapprox(result.initial.total, result.final.total; rtol=1e-5, atol=0.0)
    end

end
