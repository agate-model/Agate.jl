using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans: Clock
using Oceananigans.Units
using Oceananigans.Fields: FunctionField
using Agate.Constructors

const year = years = 365day

@testset "mass_balance" begin
    @testset "size_structured_NPZD box model" begin
        model = construct_size_structured_NPZD()

        bgc_model = Biogeochemistry(
            model(); light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
        )
        box_model = BoxModel(; biogeochemistry=bgc_model)
        set!(box_model; N=7, P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, D=0.0)

        function estimate_mass(box_model)
            return box_model.fields.N.data[1, 1, 1] +
                   box_model.fields.P1.data[1, 1, 1] +
                   box_model.fields.P2.data[1, 1, 1] +
                   box_model.fields.Z1.data[1, 1, 1] +
                   box_model.fields.Z2.data[1, 1, 1] +
                   box_model.fields.D.data[1, 1, 1]
        end

        initial_nitrogen_mass = estimate_mass(box_model)

        for n in 1:1000
            time_step!(box_model, 1.0)
        end

        final_nitrogen_mass = estimate_mass(box_model)
        @test isapprox(initial_nitrogen_mass, final_nitrogen_mass)
    end
    @testset "thunder_egg_1 box model" begin
        model = construct_thunder_egg_1()

        bgc_model = Biogeochemistry(
            model(); light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
        )
        box_model = BoxModel(; biogeochemistry=bgc_model)
        set!(box_model; DIN=7, PO4=3, P1=0.01, P2=0.01, Z1=0.05, Z2=0.05, D=0.0)

        function estimate_carbon_mass(box_model)
            return box_model.fields.DIC.data[1, 1, 1] +
                   box_model.fields.P1.data[1, 1, 1] +
                   box_model.fields.P2.data[1, 1, 1] +
                   box_model.fields.Z1.data[1, 1, 1] +
                   box_model.fields.Z2.data[1, 1, 1] +
                   box_model.fields.POC.data[1, 1, 1] +
                   box_model.fields.DOC.data[1, 1, 1]
        end

        initial_carbon_mass = estimate_carbon_mass(box_model)

        for n in 1:1000
            time_step!(box_model, 1.0)
        end

        final_carbon_mass = estimate_carbon_mass(box_model)
        @test isapprox(initial_carbon_mass, final_carbon_mass)
    end
end
