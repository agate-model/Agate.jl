using Agate
using Agate.Constructor: construct
using Agate.Library.Light
using Agate.Models: NiPiZDFactory, DarwinFactory
using Agate.Utils: box_model_mass_balance

using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans

@testset "mass_balance" begin

    function build_box_model(bgc_instance)
        bgc_model = Biogeochemistry(
            bgc_instance; light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
        )
        return BoxModel(; biogeochemistry=bgc_model)
    end

    @testset "NiPiZD box model" begin
        bgc_instance = construct(NiPiZDFactory())
        box_model = build_box_model(bgc_instance)
        set!(box_model; N=7, P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, D=0.0)

        budgets = (
            total = [:N => 1, :P1 => 1, :P2 => 1, :Z1 => 1, :Z2 => 1, :D => 1],
        )

        result = box_model_mass_balance(box_model, budgets; dt=0.1, nsteps=1000)
        @test isapprox(result.initial.total, result.final.total; rtol=1e-12, atol=0.0)
    end

    @testset "DARWIN model" begin
        bgc_instance = construct(DarwinFactory())
        box_model = build_box_model(bgc_instance)
        set!(
            box_model;
            DIN=7,
            PO4=3,
            P1=0.01,
            P2=0.01,
            Z1=0.05,
            Z2=0.05,
            DOC=0.0,
            POC=0.0,
            DON=0.0,
            PON=0.0,
            DOP=0.0,
            POP=0.0,
        )

        n2c = bgc_instance.parameters.nitrogen_to_carbon
        p2c = bgc_instance.parameters.phosphorus_to_carbon

        budgets = (
            carbon = [:DIC => 1, :P1 => 1, :P2 => 1, :Z1 => 1, :Z2 => 1, :POC => 1, :DOC => 1],
            nitrogen = [
                :DIN => 1,
                :P1 => n2c,
                :P2 => n2c,
                :Z1 => n2c,
                :Z2 => n2c,
                :PON => 1,
                :DON => 1,
            ],
            phosphorus = [
                :PO4 => 1,
                :P1 => p2c,
                :P2 => p2c,
                :Z1 => p2c,
                :Z2 => p2c,
                :POP => 1,
                :DOP => 1,
            ],
        )

        result = box_model_mass_balance(box_model, budgets; dt=1.0, nsteps=1000)
        rtol = 1e-6
        @test isapprox(result.initial.carbon, result.final.carbon; rtol=rtol)
        @test isapprox(result.initial.nitrogen, result.final.nitrogen; rtol=rtol)
        @test isapprox(result.initial.phosphorus, result.final.phosphorus; rtol=rtol)
    end
end
