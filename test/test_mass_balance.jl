using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans.Units

const year = years = 365day

@testset "mass_balance" begin
    @testset "NiPiZD box model" begin
        bgc_type = Agate.Models.NiPiZD.construct()
        model = Agate.Models.NiPiZD.instantiate(bgc_type)

        bgc_model = Biogeochemistry(
            model; light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
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

        initial_mass = estimate_mass(box_model)

        for _ in 1:1000
            time_step!(box_model, 0.1)
        end

        final_mass = estimate_mass(box_model)
        @test isapprox(initial_mass, final_mass)
    end

    @testset "DARWIN model" begin
        if isdefined(Agate.Models, :DARWIN)
            model = Agate.Models.DARWIN.construct()
            bgc_instance = model()

            bgc_model = Biogeochemistry(
                bgc_instance; light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
            )
            box_model = BoxModel(; biogeochemistry=bgc_model)
            set!(
                box_model;
                DIN=7,
                PO4=3,
                P1=0.01,
                P2=0.01,
                Z1=0.05,
                Z2=0.05,
                DOC=0,
                POC=0.0,
                DON=0,
                PON=0.0,
                DOP=0,
                POP=0.0,
            )

            function estimate_carbon_mass(box_model)
                return box_model.fields.DIC.data[1, 1, 1] +
                       box_model.fields.P1.data[1, 1, 1] +
                       box_model.fields.P2.data[1, 1, 1] +
                       box_model.fields.Z1.data[1, 1, 1] +
                       box_model.fields.Z2.data[1, 1, 1] +
                       box_model.fields.POC.data[1, 1, 1] +
                       box_model.fields.DOC.data[1, 1, 1]
            end

            function estimate_nitrogen_mass(box_model)
                return box_model.fields.DIN.data[1, 1, 1] +
                       box_model.fields.P1.data[1, 1, 1] *
                       bgc_instance.parameters.nitrogen_to_carbon +
                       box_model.fields.P2.data[1, 1, 1] *
                       bgc_instance.parameters.nitrogen_to_carbon +
                       box_model.fields.Z1.data[1, 1, 1] *
                       bgc_instance.parameters.nitrogen_to_carbon +
                       box_model.fields.Z2.data[1, 1, 1] *
                       bgc_instance.parameters.nitrogen_to_carbon +
                       box_model.fields.PON.data[1, 1, 1] +
                       box_model.fields.DON.data[1, 1, 1]
            end

            function estimate_phosphorus_mass(box_model)
                return box_model.fields.PO4.data[1, 1, 1] +
                       box_model.fields.P1.data[1, 1, 1] *
                       bgc_instance.parameters.phosphorus_to_carbon +
                       box_model.fields.P2.data[1, 1, 1] *
                       bgc_instance.parameters.phosphorus_to_carbon +
                       box_model.fields.Z1.data[1, 1, 1] *
                       bgc_instance.parameters.phosphorus_to_carbon +
                       box_model.fields.Z2.data[1, 1, 1] *
                       bgc_instance.parameters.phosphorus_to_carbon +
                       box_model.fields.POP.data[1, 1, 1] +
                       box_model.fields.DOP.data[1, 1, 1]
            end

            initial_carbon = estimate_carbon_mass(box_model)
            initial_nitrogen = estimate_nitrogen_mass(box_model)
            initial_phosphorus = estimate_phosphorus_mass(box_model)

            for _ in 1:1000
                time_step!(box_model, 1)
            end

            final_carbon = estimate_carbon_mass(box_model)
            final_nitrogen = estimate_nitrogen_mass(box_model)
            final_phosphorus = estimate_phosphorus_mass(box_model)

            rtol = 1e-6
            @test isapprox(initial_carbon, final_carbon, rtol=rtol)
            @test isapprox(initial_nitrogen, final_nitrogen, rtol=rtol)
            @test isapprox(initial_phosphorus, final_phosphorus, rtol=rtol)
        else
            @test true
        end
    end
end
