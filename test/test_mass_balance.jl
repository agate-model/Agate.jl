using Agate
using Agate.Library.Light
using OceanBioME
using OceanBioME: Biogeochemistry
using Oceananigans
using Oceananigans: Clock
using Oceananigans.Units
using Oceananigans.Fields: FunctionField
using Agate.Models

const year = years = 365day

@testset "mass_balance" begin
    @testset "size_structured_NPZD box model" begin
        N2P2ZD_constructed = Agate.Models.NiPiZD.construct()
        model = Agate.Models.NiPiZD.instantiate(N2P2ZD_constructed)

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

        initial_nitrogen_mass = estimate_mass(box_model)

        for n in 1:1000
            time_step!(box_model, 0.1)
        end

        final_nitrogen_mass = estimate_mass(box_model)
        @test isapprox(initial_nitrogen_mass, final_nitrogen_mass)
    end

    @testset "Agate.jl-DARWIN model" begin
        model = Agate.Models.DARWIN.construct()

        bgc_model = Biogeochemistry(
            model(); light_attenuation=FunctionFieldPAR(; grid=BoxModelGrid())
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
                   box_model.fields.P1.data[1, 1, 1] * model().nitrogen_to_carbon +
                   box_model.fields.P2.data[1, 1, 1] * model().nitrogen_to_carbon +
                   box_model.fields.Z1.data[1, 1, 1] * model().nitrogen_to_carbon +
                   box_model.fields.Z2.data[1, 1, 1] * model().nitrogen_to_carbon +
                   box_model.fields.PON.data[1, 1, 1] +
                   box_model.fields.DON.data[1, 1, 1]
        end

        function estimate_phosphorus_mass(box_model)
            return box_model.fields.PO4.data[1, 1, 1] +
                   box_model.fields.P1.data[1, 1, 1] * model().phosphorus_to_carbon +
                   box_model.fields.P2.data[1, 1, 1] * model().phosphorus_to_carbon +
                   box_model.fields.Z1.data[1, 1, 1] * model().phosphorus_to_carbon +
                   box_model.fields.Z2.data[1, 1, 1] * model().phosphorus_to_carbon +
                   box_model.fields.POP.data[1, 1, 1] +
                   box_model.fields.DOP.data[1, 1, 1]
        end

        carbon_mass_records = Float64[]
        nitrogen_mass_records = Float64[]
        phosphorus_mass_records = Float64[]
        initial_carbon_mass = estimate_carbon_mass(box_model)
        initial_nitrogen_mass = estimate_nitrogen_mass(box_model)
        initial_phosphorus_mass = estimate_phosphorus_mass(box_model)
        push!(carbon_mass_records, initial_carbon_mass)
        push!(nitrogen_mass_records, initial_nitrogen_mass)
        push!(phosphorus_mass_records, initial_phosphorus_mass)

        for n in 1:1000
            time_step!(box_model, 1)
        end

        final_carbon_mass = estimate_carbon_mass(box_model)
        final_nitrogen_mass = estimate_nitrogen_mass(box_model)
        final_phosphorus_mass = estimate_phosphorus_mass(box_model)

        rtol = 1e-6

        @test isapprox(initial_carbon_mass, final_carbon_mass, rtol=rtol)
        @test isapprox(initial_phosphorus_mass, final_phosphorus_mass, rtol=rtol)
        @test isapprox(initial_nitrogen_mass, final_nitrogen_mass, rtol=rtol)
    end
end
