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

    @testset "Generic multi-nutrient conservation" begin
        grid = BoxModelGrid()
        bgc_instance = multi_nutrient_test_model()
        box_model = build_box_model(bgc_instance, grid)
        set!(
            box_model;
            DIC=2.0,
            DIN=7.0,
            PO4=3.0,
            DOC=0.0,
            POC=0.0,
            DON=0.0,
            PON=0.0,
            DOP=0.0,
            POP=0.0,
            P_1=0.01,
        )

        n2c = bgc_instance.parameters.nitrogen_to_carbon
        p2c = bgc_instance.parameters.phosphorus_to_carbon
        budgets = (
            carbon=[:DIC => 1, :P_1 => 1, :POC => 1, :DOC => 1],
            nitrogen=[:DIN => 1, :P_1 => n2c, :PON => 1, :DON => 1],
            phosphorus=[:PO4 => 1, :P_1 => p2c, :POP => 1, :DOP => 1],
        )

        result = box_model_mass_balance(box_model, budgets; dt, nsteps)
        @test isapprox(result.initial.carbon, result.final.carbon; rtol=1e-6)
        @test isapprox(result.initial.nitrogen, result.final.nitrogen; rtol=1e-6)
        @test isapprox(result.initial.phosphorus, result.final.phosphorus; rtol=1e-6)
    end

end
